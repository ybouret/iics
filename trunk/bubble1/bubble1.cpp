#include "arrays.h"
#include <mpi.h>
/*** For Silo ***/
#include <silo.h>
#include <stdio.h> 
/* Instrumenting with visit */
 
#define L 50.0
#define ANX 100
#define pi 3.141569
static const        size_t NC_s = 2;     /*!< number of scalar components */
static const        size_t NC_v = 1;     /*!< number  vector component */
static const        size_t NC_t = 1;     /*!< number of temp scalar component for fieldst */
size_t              NC;                  /* Total number of component NC=Nc_S+2*NC_v*/
const char          *cpntName_s[] = {"pressure","c"};  /*Name of the components*/
const char          *cpntName_v[] = {"velocity"};  /*Name of the components*/
static indx_t       Nx = ANX;   /*!< 0 -> Nx-1      */
static indx_t       Ny = 1;   /*!< 0 -> Ny-1      */
static indx_t       Nz = ANX;   /*!< 0 -> Nz-1      */
static indx_t       NG = 1;     /*!< #ghosts        */
static real_t       Lx = L;
static real_t       Ly = L;
static real_t       Lz = L; 

static real_t   ****fields          = NULL;
static real_t   ****fieldst         = NULL; // temporary scalar field
static real_t   ****dfields         = NULL;
static int          size            = 0; /*!< #procs     */
static int          rank            = 0; /*!< proc id    */
static int          above           = 0; /*!< proc above */
static int          below           = 0; /*!< proc below */

static indx_t       xmin            = 0;
static indx_t       xmax            = 0;
static indx_t       ymin            = 0;
static indx_t       ymax            = 0;
static size_t       items_per_slice = 0;
static size_t       items_per_field = 0;
static indx_t       zmin            = 0;/*!< without ghosts  */  
static indx_t       zmax            = 0; /*!< without ghosts */
static indx_t       zlo             = 0; /*!< with ghosts    */
static indx_t       zhi             = 0; /*!< with ghosts    */

real_t       dx=1,dy=1,dz=1;
real_t       idx=1.0/dx,idy=1.0/dy,idz=1.0/dz;
static real_t       idx2=idx*idx,idy2=idy*idy,idz2=idz*idz;
static real_t       dt = 0.01;
float               *rmesh[3];


#define _CHECK(MPI_ROUTINE) do { const int __rc = MPI_ROUTINE; if( MPI_SUCCESS != __rc ) { fprintf( stderr, "Failure: " #MPI_ROUTINE "\n" ); exit(-1); } } while(0)
#define _BARRIER _CHECK(MPI_Barrier(MPI_COMM_WORLD))
static MPI_Request *requests = NULL;
static size_t       num_reqs = 0;
static const int    diff_tag = 7;

#define alea (rand()/(0.0+RAND_MAX)*2-1)

int   rmesh_dims[3];



#define eta 0.05
#define gamma 0.1
#define noise 0.1
#define NB 200
typedef struct { 
    real_t x[NB];
    real_t z[NB]; 
    real_t pressure;
    int domain[NB];
    int  n; 
}          Bubble; 
MPI_Datatype bubbletype, oldtypes[2];  
Bubble b0;
Bubble *bt;


#include "stubs.c"
#include "ui.cpp"
#include "writeToSilo.c"
#include "mpiPartition.c"

void initiateBubble()
{
    MPI_Status stat; 
    /* MPI_Aint type used to be consistent with syntax of */ 
    /* MPI_Type_extent routine */ 
    MPI_Aint    offsets[2], extent; 
    int          blockcounts[2]; 

    
    /* Setup description of the 4 MPI_FLOAT fields x, y, z, velocity */ 
    offsets[0] = 0; 
    oldtypes[0] = MPI_DOUBLE; 
    blockcounts[0] = 1+2*NB; 
    MPI_Type_extent(MPI_DOUBLE, &extent);
    
    offsets[1] = blockcounts[0] * extent; 
    oldtypes[1] = MPI_INT; 
    blockcounts[1] = 1+NB; 
    
    /* Now define structured type and commit it */ 
    MPI_Type_struct(2, blockcounts, offsets, oldtypes, &bubbletype); 
    MPI_Type_commit(&bubbletype);
    
    
    //if(rank==0)
    bt=(Bubble *) malloc(size*sizeof(Bubble));
    
    /* Initialize the bubble and then send it to each core */ 
    b0.pressure=-1;
    if(rank==0)
    {
        indx_t i;
        //radius and pressure of the bubble
        real_t R=10;
        real_t Pg=5;
        
        for(i=0;i<NB;i++)
        {
            b0.x[i]=R*cos(2*3.1415*i/(NB-1.0))+Lx*0.5;
         //   b0.x[i]=10;
            b0.z[i]=R*sin(2*3.1415*i/(NB-1.0))+Lz*0.5;
           // b0.z[i]=i/(NB-1.0)*Lz;
          //  b0.z[i]=i/(Nx-1.0)*Lz;
            b0.domain[i]=-1;
        
        }
        b0.pressure=Pg;
        b0.n=NB;
    }
    MPI_Bcast(&b0,1,bubbletype,0,MPI_COMM_WORLD);
    /*
    if(rank==0)
    {
        int i;

        for(i=0;i<Nx;i++)
            fprintf(stderr,"\n z[%d]=%g is in domain %d\n",i,b0.z[i], (int)floor(b0.z[i]/Lz*size));
    }
  */
   
}
double f(real_t x)
{
    return x-x*x*x;
}
void putBoundaryOnPressure()
{
    int i,k;
    real_t ***pt=fieldst[0];
    real_t ***p=fields[0];
    real_t ***c=fields[1];
    
    for(i=xmin+1;i<=xmax-1;i++)
    {
        for(k=zmin;k<=zmax;k++)  
        {
            if(pt[k][0][i]>0)
            {
                p[k][0][i]=pt[k][0][i];
                c[k][0][i]=pt[k][0][i];

            }
            
        }
    }
    
}
double gaussSeidel(int color,real_t omega)
{
	indx_t i,ip,im,j=0,k,km,kp;
    real_t ***p=fields[0];
    real_t ***pt=fieldst[0];
    real_t anorm=0.0,resid=0.0;
    real_t x,z;
    real_t R=5;
    real_t Pg=0;
    real_t a=1.0/(dz*dz);
    real_t b=1.0/(dx*dx);
    real_t c=-2.0*(a+b);
    real_t ic=1.0/c;
    for(k=zmin;k<=zmax;k++)  
    {
       // p[k][j][xmax]=p[k][j][xmax-1];
        p[k][j][xmax]=1;
       // p[k][j][xmin]=p[k][j][xmin+1];
        p[k][j][xmin]=0;
    }
    
    putBoundaryOnPressure();
    
    for(i=xmin+1;i<=xmax-1;i++)
    {
        ip=i+1;
        im=i-1;
        x=i*dx-0.5*Lx;

        for(k=zmin;k<=zmax;k++)  
        {
            km=k-1;
            kp=k+1;
            z=k*dz-0.5*Lz;

            if((i+k)%2==color)
            {                
                if(pt[k][j][i]>0)
                {
                    p[k][j][i]=pt[k][j][i];
                    resid=0;
                }
                else
                {

                    resid=( a*(p[kp][j][i]+p[km][j][i])
                                 +b*(p[k][j][im]+p[k][j][ip])
                                 +c*p[k][j][i]);
                    p[k][j][i]-=omega*resid*ic;     
                }
                anorm+=fabs(resid);

            }
        }
    }
    
    sendRequests(0);
    waitRequests(0);
    return anorm;
}

void computePressure(int nbIter,real_t eps)
{
    real_t myError=0;
    real_t rho_s=(cos(pi/(Nx-1))+(dx/dz)*(dx/dz)*(pi/(Nz-1)))/(1+(dx/dz)*(dx/dz));
    real_t rho_s2=rho_s*rho_s;

    //Chebychev acceleration
    real_t omega=1.0;

    double error=1.0;
    indx_t k,i;
    

    for(k=0;(k<nbIter)&& (error>eps);k++)
    {
        
        myError=gaussSeidel(0,omega);
        if(k==0)
            omega=1.0/(1-rho_s2*0.5);
        else
            omega=1.0/(1-rho_s2*0.25*omega);
        myError+=gaussSeidel(1,omega);
        omega=1.0/(1-rho_s2*0.25*omega);
        MPI_Reduce(&myError, &error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
        MPI_Bcast( &error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    }
    if(rank==0)
    {
        if(error>eps)
            fprintf(stderr,"\n not converged after %d iterations\t  total Error is %g\n",(int) k,error);
    }
}

double gaussNLSeidel(int color)
{
	indx_t i,ip,im,j=0,k,km,kp,part;
    real_t ***p=fields[0];
    real_t ***pt=fieldst[0];
    real_t anorm=0.0,resid=0.0;
    real_t x,z;
    real_t R=10;
    real_t Pg=5;
    real_t a=1.0/(dz*dz);
    real_t b=1.0/(dx*dx);
    real_t c=-2.0*(a+b);
    real_t ic=1.0/c;
    real_t ap,epsilon=0.0001;
    //    fprintf(stderr,"proc %d: xmin=%d\t xmax=%d\n",rank,(int)xmin,(int)xmax);
    for(k=zmin;k<=zmax;k++)  
    {
        p[k][j][xmax]=1;
        p[k][j][xmin]=p[k][j][xmin+1];
       // p[k][j][xmin]=0;
    }
   // putBoundaryOnPressure();
    
    for(i=xmin+1;i<=xmax-1;i++)
    {
        ip=i+1;
        im=i-1;
        x=i*dx-0.5*Lx;
        
        for(k=zmin;k<=zmax;k++)  
        {
            km=k-1;
            kp=k+1;
            z=k*dz-0.5*Lz;
            
            if((i+k)%2==color)
            {                
                if(pt[k][0][i]>0)
                {
                    p[k][j][i]=pt[k][0][i];
                    resid=0;
                }
                else
                {
                    ap=p[k][j][i];

                    
                    resid=( a*(p[kp][j][i]+p[km][j][i])
                           +b*(p[k][j][im]+p[k][j][ip])
                           +c*p[k][j][i]+f(ap));
                    p[k][j][i]-=resid/(c+(f(ap+epsilon)-f(ap-epsilon))/(2*epsilon));     
                }
                anorm+=fabs(resid);
                
            }
        }
    }
    
    sendRequests(0);
    waitRequests(0);
    return anorm;
}
void computeNLPressure(int nbIter,real_t eps)
{
    real_t myError=0;
    real_t rho_s=(cos(pi/(Nx-1))+(dx/dy)*(dx/dy)*(pi/(Nz-1)))/(1+(dx/dy)*(dx/dy));
    real_t rho_s2=rho_s*rho_s;
    
    //Chebychev acceleration
    real_t omega=1.0;
    
    double error=1.0;
    indx_t k,i;
    
    
    for(k=0;(k<nbIter)&& (error>eps);k++)
    {
        
        myError=gaussNLSeidel(0);
        myError+=gaussNLSeidel(1);
        MPI_Reduce(&myError, &error, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
        MPI_Bcast( &error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
        if(rank==0)
            fprintf(stderr,"\n not converged after %d iterations\t  total Error is %g\n",(int) k,error);
        
    }
    if(rank==0)
    {
        if(error>eps)
            fprintf(stderr,"\n not converged after %d iterations\t  total Error is %g\n",(int) k,error);
            }
}
void computeVelocity()
{
    indx_t k,km,kp,i,ip,im,j=0;
    
    real_t ***p=fields[0];
    real_t ***u=fields[2];
    real_t ***v=fields[3];

    for(k=zmin;k<=zmax;k++)  
    {
        km=k-1;
        kp=k+1;
        for(i=xmin;i<=xmax;i++)
        {
            ip=i+1;
            im=i-1;
            if(i==0)
            {
                u[k][j][i]=0;
                v[k][j][i]=-(p[kp][j][i]-p[km][j][i])/dz*0.5;
            }
            else if(i==Nx-1)
            {
                u[k][j][i]=-(1-p[k][j][im])/dx*0.5;
                v[k][j][i]=-(p[kp][j][i]-p[km][j][i])/dz*0.5;

            }
            else
            {
                u[k][j][i]=-(p[k][j][ip]-p[k][j][i])/dx*0.5;
                v[k][j][i]=-(p[kp][j][i]-p[k][j][i])/dz*0.5;
            } 
           // u[k][j][i]=1.0;
           // v[k][j][i]=0.1;
        }
    }
    sendRequests(2);
    sendRequests(3);
    waitRequests(2);
    waitRequests(3);

}
void advectBubble(simulation_data *sim)
{
    
    int i,k;
    int ip,kp;
    int ipp,ipm;
   // int p;
    int particleDomain;
    real_t ***pt=fieldst[0];
    real_t ***p=fields[0];
    real_t ***c=fields[1];
    real_t ***u=fields[2];
    real_t ***v=fields[3];
    real_t cosTheta;
    real_t sinTheta;
    real_t idis, tx,tz;
    MPI_Status stat; 
    
    for(i=xmin;i<=xmax;i++)
        for(k=zmin-1;k<=zmax+1;k++)
        {
            pt[k][0][i]=0;
            c[k][0][i]=0;
        }
    
    dt=0.05;

    for(i=0;i<NB;i++)
    {
        //b0.domain[i]=-2;
        particleDomain=(int)floor(b0.z[i]/Lz*size);
        
        if(particleDomain==rank)
        {
            indx_t kshift;
            indx_t ishift;
            
            b0.domain[i]=particleDomain;

            ip=(int)(b0.x[i]/Lx*(Nx));
            if(ip<0 ||(ip>xmax))
                    fprintf(stderr,"sortie de tableau x dans advecBubble ip=%d\n", (int) ip);

            kp=(int)(b0.z[i]/Lz*(Nz));
            if(kp<zmin-1 ||(kp>zmax+1))
                fprintf(stderr,"sortie de tableau z dans advecBubble : %d<=kp=%d<=%d\n",(int) zmin,(int) kp,(int) zmax);
            
            // interpolated velocity
            b0.x[i]+=dt*(u[kp][0][ip-1]+(b0.x[i]-(ip-1)*dx)/dx*(u[kp][0][ip]-u[kp][0][ip-1]));
            b0.z[i]+=dt*(v[kp-1][0][ip]+(b0.z[i]-(kp-1)*dz)/dz*(v[kp][0][ip]-v[kp][0][ip]));
            
            // we impose the pressure inside the bubble in order to avoid numerical problems.
            ipp=i+1;
            ipm=i-1;
            if(i==0)
                ipm=NB-1;
            if(i==NB-1)
                ipp=0;
        
            tx=(b0.x[ipp]-b0.x[ipm]);
            tz=(b0.z[ipp]-b0.z[ipm]);
            idis=sqrt(tx*tx+tz*tz);
            if(idis==0)
                fprintf(stderr,"error in advectBubble: two particles are at the same place\n");
            else
                idis=1.0/idis;
            cosTheta=tx*idis;
            kshift=
            sinTheta=tz*idis;
        //    fprintf(stderr,"particule %d cos=%d,sin=%d\n",i,(int) lrint(cosTheta),(int) lrint(sinTheta));
            //pt[kp-2*((int) cosTheta)][0][ip+2*((int) -sinTheta)]=b0.pressure;
            //c[kp-((int) lrint(cosTheta))][0][ip-((int)lrint(-sinTheta))]=b0.pressure;
            //c[kp-((int) lrint(cosTheta))][0][ip-((int)lrint(-sinTheta))]=b0.pressure;
            kshift=((int) lrint(cosTheta));
            ishift=((int)lrint(-sinTheta));
                if((zmin<=kp+kshift)&&(kp+kshift<=zmax))
            {
                pt[kp+kshift][0][ip+ishift]=b0.pressure;
                c[kp+kshift][0][ip+ishift]=b0.pressure;

            }
            else
            {
                pt[kp][0][ip+ishift]=b0.pressure;
                c[kp][0][ip+ishift]=b0.pressure;
            }

            
        }
        else
        {
            b0.domain[i]=-1;
        }
    }
    /***** problem avec le changement de domaine !!!!!*/
   _CHECK(MPI_Gather(&b0, 1, bubbletype, bt, 1, bubbletype, 0,MPI_COMM_WORLD));
    if(rank==0)
    {
        for(i=0;i<NB;i++)
        {
            for(k=0;k<size;k++)
            {  
               // if(b0.domain[i]==-2)
                 //   exit(0);
                if(bt[k].domain[i]>=0)
                {
                    
                   // particleDomain=(int)floor(b0.z[i]/Lz*size);
                    b0.domain[i]=bt[k].domain[i];
                    b0.x[i]=bt[k].x[i];
                    b0.z[i]=bt[k].z[i];                     
                }
            }
           
        }
    }
    
    _CHECK(MPI_Bcast(&b0,1,bubbletype,0,MPI_COMM_WORLD));
}
static void integrate(simulation_data *sim)
{
    indx_t i,j;
    
    computeVelocity();
    advectBubble(sim);
    
    computePressure(100000,1e-4);
    
 //   computeNLPressure(100,1e-3);
    
    /*
    
    //euler integration
	for( i=0; i < NC; ++i )
	{
		real_t ***f = fields[i];
		real_t ***df = dfields[i];
        
        real_t       *dst = &f[zmin][ymin][xmin];
        const real_t *src = &df[zmin][ymin][xmin];
        for( j=0; j < items_per_field; ++j )
        {
            dst[j] += dt * (src[j]);
        }
    }
     */
}


static void init_fields(int cond)
{
	indx_t i,j=0,k;
    real_t x,y,z;
    real_t ***p=fields[0];

    srand(time(NULL)+rank);
	
    if(rank==0)
        fprintf(stderr,"initial condition is %d\n",cond);
    for(k=zmax;k>=zmin;--k)  
	{
        z=k*dz-Lz*0.5;
	//	for(j=ymax;j>=ymin;--j)
		{
            y=j*dy-Ly*0.5;
			for(i=xmax;i>=xmin;--i)
			{
                x=i*dx-Lx*0.5;
              //  if(x*x+z*z<100)
                switch(cond)
                {
                    case 0:   
                        p[k][j][i] =2;
                        if(i==xmax)
                            p[k][j][i] =1.0;
                        else if(i==0)
                            p[k][j][i] =0.0;

                        break;
                
                    default:   
                        p[k][j][i] =0.6+0.000*alea;
                        break;
                }
               
            }
		}
	}
    sendRequests(0);
    waitRequests(0);
}


void simulate_one_timestep(simulation_data *sim)
{
    /* simulate 1 time step. */
    int i;
    double elapsedTime;
    ++sim->cycle;
    sim->time += dt;
    if(rank==0)
        elapsedTime=MPI_Wtime();
    
    // Diffusion 
    for(i=0;i<1;i++)
    {
       integrate(sim);
    }
    
    if(sim->savingFiles==1)
    {
        writeDomain2D(sim->cycle);
        write_master(sim->cycle);
    }
    
    
    if(sim->visitIsConnected==1)
    {
        VisItTimeStepChanged();
        VisItUpdatePlots();
    }
    
    if(rank==0)  
    {
        fprintf(stderr,"%d cores:simulating: cycle=%d, time=%lg in %g s\n",sim->par_size, sim->cycle, sim->time,MPI_Wtime()-elapsedTime);
        fflush(stderr);
    }
}
int main(int argc, char *argv[] )
{
	//unsigned count = 0;
    //int toWrite=0;
    //double startTime,endTime;
    
    
    
	/***************************************************************************
	 * geometry setup
	 **************************************************************************/
	xmin            = 0;
	xmax            = Nx-1;
	ymin            = 0;
	ymax            = Ny-1;
	items_per_slice = Nx * Ny;
	
	dx   = Lx / ( Nx-1.0 );
	dy   = Ly / ( Ny-1.0 );
	dz   = Lz / ( Nz-1.0 );
	idx2 = 1/(dx*dx); 
	idy2 = 1/(dy*dy);
	idz2 = 1/(dz*dz);
	srand( time(NULL) );
	setvbuf( stderr, NULL, 0, _IONBF );
    
    
	
	_CHECK(MPI_Init(&argc,&argv));
	/***************************************************************************
	 * init of MPI, slicing of the domains, etc....
	 **************************************************************************/
    NC=NC_s+2*NC_v;
    initSimulation(); 
    
    
    /***************************************************************************
     * Initialize environment variables. 
	 **************************************************************************/
    

    
    
    simulation_data sim;
    simulation_data_ctor(&sim,rank,size);
    //  VisItSetDirectory("/Applications/VisIt.app/Contents/Resources/2.4.0/darwin-x86_64");
    VisItSetupEnvironment();
    
    /***************************************************************************
     * Install callback functions for global communication.. 
	 **************************************************************************/
    VisItSetBroadcastIntFunction(visit_broadcast_int_callback);
    VisItSetBroadcastStringFunction(visit_broadcast_string_callback);
    
    /***************************************************************************
     * Tell visit whether the simulation is parallel. 
	 **************************************************************************/
    VisItSetParallel(sim.par_size > 1);
    VisItSetParallelRank(sim.par_rank);   
    if(rank==0)
    {
        
        if(VisItInitializeSocketAndDumpSimFile("bubble1",
                                               "Parallel C prototype simulation connects to VisIt",
                                               argv[0], NULL, "bubble1.ui", NULL))
      
        {
            fprintf(stderr,"VisItInitializeSocketAndDumpSimFile ok\n");
        }
        else
        {
            fprintf(stderr,"problem with VisItInitializeSocketAndDumpSimFile \n");
        }
    }
	init_fields(0);
    initiateBubble();
    sim.savingFiles=1;
    if(rank==0) VisItOpenTraceFile("./TraceFileOfLibSim.txt");
    set_interface(&sim);

    mainloop(&sim);
    if(rank==0) VisItCloseTraceFile();
    /***************************************************************************
	 * measuring the communication time
	 **************************************************************************/
    //fprintf(stderr,"Communication time of rank %d is %g\n",rank,mesureTimeForExhangingGhost);
    
	/***************************************************************************
	 * simulation
	 **************************************************************************
	 _BARRIER;
     startTime=MPI_Wtime();
     for( count=0; count < 100; ++count )
     {
     simulate_one_timestep(&sim);
     }
     _BARRIER;
     endTime=MPI_Wtime();
     if(rank==0)
     fprintf(stderr,"rank=%d\t ran %d counts in %g\n",rank,count,endTime-startTime);
     ***************************************************************************
	 * MPI cleanup
	 **************************************************************************/
	simulation_data_dtor(&sim);
    MPI_Finalize();
    
	return 0;
}
