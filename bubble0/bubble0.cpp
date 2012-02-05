#include "arrays.h"
#include <mpi.h>
/*** For Silo ***/
#include <silo.h>
#include <stdio.h> 
/* Instrumenting with visit */




 

static const size_t NC = 4;     /*!< two components */
const char                *cpntName[] = {"rho","u","v","pressure"};  /*Name of the components*/
static indx_t       Nx = 128;   /*!< 0 -> Nx-1      */
static indx_t       Ny = 1;   /*!< 0 -> Ny-1      */
static indx_t       Nz = 128;   /*!< 0 -> Nz-1      */
static indx_t       NG = 2;     /*!< #ghosts        */
static real_t       Lx = 128.0;
static real_t       Ly = 64.0;
static real_t       Lz = 128.0; 

static real_t   ****fields          = NULL;
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

static real_t       dx=1,dy=1,dz=1;
static real_t       idx=1.0/dx,idy=1.0/dy,idz=1.0/dz;
static real_t       idx2=idx*idx,idy2=idy*idy,idz2=idz*idz;
static real_t       dt = 0.01;

#define _CHECK(MPI_ROUTINE) do { const int __rc = MPI_ROUTINE; if( MPI_SUCCESS != __rc ) { fprintf( stderr, "Failure: " #MPI_ROUTINE "\n" ); exit(-1); } } while(0)
#define _BARRIER _CHECK(MPI_Barrier(MPI_COMM_WORLD))
static MPI_Request *requests = NULL;
static size_t       num_reqs = 0;
static const int    diff_tag = 7;

#define alea (rand()/(0.0+RAND_MAX)*2-1)

int   rmesh_dims[3];
static void initMetaBubble();
#include "stubs.c"
#include "ui.cpp"
#include "writeToSilo.c"


static void create_fields()
{
	/***************************************************************************
	 ** NC component + 1 field for the laplacian value
	 **************************************************************************/
	size_t i;
	fields  = (real_t ****)calloc(NC,sizeof(real_t ***));
	dfields = (real_t ****)calloc(NC,sizeof(real_t ***));
	if( !fields || !dfields )
	{
		perror("memory allocation error in fields level-1");
		exit(-1);
	}
    
	for(i=0;i<NC;++i)
	{
		dfields[i] = icp_create_array3D(xmin,xmax,ymin,ymax,zlo,zhi);
		fields[i]  = icp_create_array3D(xmin,xmax,ymin,ymax,zlo,zhi);
		if( !fields[i] || !dfields[i])
		{
			while( i > 0 )
			{
				icp_delete_array3D(fields[--i],xmin,xmax,ymin,ymax,zlo,zhi);
				icp_delete_array3D(dfields[--i],xmin,xmax,ymin,ymax,zlo,zhi);
			}
			free(dfields);
			dfields=NULL;
			perror("fields level-2");
			exit(-1);
		}
	}
}

static void delete_fields()
{
	if( fields)
	{
		int i=NC;
		while( i > 0 )
		{
			icp_delete_array3D(fields[--i],xmin,xmax,ymin,ymax,zlo,zhi);
			icp_delete_array3D(dfields[--i],xmin,xmax,ymin,ymax,zlo,zhi);
		}
		free(fields);
		free(dfields);
		fields = NULL;
		dfields = NULL;
	}
}


static void create_requests()
{
	//const size_t nitems = items_per_slice * NG;
	//size_t       i;
	
	num_reqs = 4 * NC;
	requests = (MPI_Request *)calloc(num_reqs,sizeof(MPI_Request));
	if( !requests )
	{
		perror("requests");
		exit(-1);
	}
	
    /*
     for( i=0; i < NC; ++i )
     {
     const size_t j = i * 4;
     // send information to below 
     _CHECK(MPI_Send_init( &fields[i][zmin][ymin][xmin],    nitems, ICP_REAL, below, diff_tag, MPI_COMM_WORLD, &requests[j+0] ));
     
     // send information to above 
     _CHECK(MPI_Send_init( &fields[i][zmax-NG][ymin][xmin], nitems, ICP_REAL, above, diff_tag, MPI_COMM_WORLD, &requests[j+1] ));
     
     // recv information from below
     _CHECK(MPI_Recv_init( &fields[i][zlo][ymin][xmin],     nitems, ICP_REAL, below, diff_tag, MPI_COMM_WORLD, &requests[j+2] ));
     
     // recv information from above 
     _CHECK(MPI_Recv_init( &fields[i][zmax+1][ymin][xmin],  nitems, ICP_REAL, above, diff_tag, MPI_COMM_WORLD, &requests[j+3] ));
     } 
     */
}

static void delete_requests()
{
	if( requests )
	{
		free( requests );
		requests = NULL;
	}
}


static void exchange_ghosts()
{
	const size_t nitems = items_per_slice * NG;
	size_t       i;
    MPI_Status status;
	
	num_reqs = 4 * NC;
	requests = (MPI_Request *)calloc(num_reqs,sizeof(MPI_Request));
	if( !requests )
	{
		perror("requests");
		exit(-1);
	}
	
	for( i=0; i < NC; ++i )
	{
		const size_t j = i * 4;
		// send information to below 
		_CHECK(MPI_Isend( &fields[i][zmin][ymin][xmin],    nitems, ICP_REAL, below, diff_tag, MPI_COMM_WORLD, &requests[j+0] ));
		
		// send information to above 
		_CHECK(MPI_Isend( &fields[i][zmax-NG][ymin][xmin], nitems, ICP_REAL, above, diff_tag, MPI_COMM_WORLD, &requests[j+1] ));
		
		// recv information from below
		_CHECK(MPI_Irecv( &fields[i][zlo][ymin][xmin],     nitems, ICP_REAL, below, diff_tag, MPI_COMM_WORLD, &requests[j+2] ));
		
		// recv information from above 
		_CHECK(MPI_Irecv( &fields[i][zmax+1][ymin][xmin],  nitems, ICP_REAL, above, diff_tag, MPI_COMM_WORLD, &requests[j+3] ));
	}
    for(i=0;i<num_reqs;i++)
        _CHECK(MPI_Wait(&requests[i],&status));
}
/*****************************************************************************************************
 *     Here we start the send/recv requests for variable i
 *****************************************************************************************************/
static void sendRequests(int i)
{
	const size_t nitems = items_per_slice * NG;
    // MPI_Status status;
    const size_t j = i * 4;
    
	if( !requests )
	{
		perror("requests in sendRequests");
		exit(-1);
	}
    
    // send information to below 
    _CHECK(MPI_Isend( &fields[i][zmin][ymin][xmin],    nitems, ICP_REAL, below, diff_tag, MPI_COMM_WORLD, &requests[j+0] ));
    
    // send information to above 
    _CHECK(MPI_Isend( &fields[i][zmax+1-NG][ymin][xmin], nitems, ICP_REAL, above, diff_tag, MPI_COMM_WORLD, &requests[j+1] ));
    
    // recv information from below
    _CHECK(MPI_Irecv( &fields[i][zlo][ymin][xmin],     nitems, ICP_REAL, below, diff_tag, MPI_COMM_WORLD, &requests[j+2] ));
    
    // recv information from above 
    _CHECK(MPI_Irecv( &fields[i][zmax+1][ymin][xmin],  nitems, ICP_REAL, above, diff_tag, MPI_COMM_WORLD, &requests[j+3] ));
	   
}
/*****************************************************************************************************
 *     Here we wait the send/recv requests for variable i
 *****************************************************************************************************/
static void waitRequests(int i)
{
    const size_t j = i * 4;
    int k;
    MPI_Status status;
    
    
    if( !requests )
	{
		perror("requests in waitRequests");
		exit(-1);
	}
    
    for(k=0;k<4;k++)
        _CHECK(MPI_Wait(&requests[j+k],&status));
}


real_t pressure(real_t x)
{
    return x-1.9*x*x+x*x*x;
//    return x;
}
#define eta 0.05
#define gamma 0.1
#define noise 0.0



real_t dxx(indx_t var,indx_t i,indx_t j,indx_t k)
{
    real_t ***f=fields[var];
    int im=i-1,ip=i+1;
    
    if( im < xmin ) im = xmax;
    if( ip > xmax ) ip = xmin;
    
    return idx2*(f[k][j][ip]-2*f[k][j][i]+f[k][j][im]);
}

real_t dzz(indx_t var,indx_t i,indx_t j,indx_t k)
{
    real_t ***f=fields[var];
    
    indx_t km=k-1; /* always valid for there are ghosts */
    indx_t kp=k+1; /* always valid for there are ghosts */

    
    return idz2*(f[kp][j][i]-2*f[k][j][i]+f[km][j][i]);
}
real_t lap(indx_t var,indx_t i,indx_t j,indx_t k)
{
    return dxx(var,i,j,k)+dzz(var,i,j,k);
}
real_t dxrho(indx_t i,indx_t j,indx_t k)
{
    real_t ***rho=fields[0];
    int im=i-1,ip=i+1;
    
    if( im < xmin ) im = xmax;
    if( ip > xmax ) ip = xmin;
    return idx*(rho[k][j][ip]-rho[k][j][im])*0.5;
}

real_t dzrho(indx_t i,indx_t j,indx_t k)
{
    real_t ***rho=fields[0];
    indx_t km=k-1; /* always valid for there are ghosts */
    indx_t kp=k+1; /* always valid for there are ghosts */
    
    return idz*(rho[kp][j][i]-rho[km][j][i])*0.5;
}
real_t sigma_xx(indx_t i,indx_t j,indx_t k)
{
    real_t ***rho =fields[0];
    real_t ***U   =fields[1];
    
    indx_t im=i-1;
    if( im < xmin ) im = xmax;
    
    return 
    -pressure(rho[k][j][i])
    +gamma*lap(0,i,j,k)
    +2.0*eta*idx*(U[k][j][i]-U[k][j][i-1]);
    //+2.0*eta*idx*(U[k][j][i]/rho[k][j][i]-U[k][j][i-1]/rho[k][j][i-1]);
}
real_t sigma_xz(indx_t i,indx_t j,indx_t k)
{
//    real_t ***rho =fields[0];
    real_t ***U   =fields[1];
    real_t ***W   =fields[2];
    
    indx_t im=i-1;
    if( im < xmin ) im = xmax;

    return eta*(
                idx*(W[k][j][i]-W[k][j][im])+
                idz*(U[k][j][i]-U[k-1][j][i])   
                //idx*(W[k][j][i]/rho[k][j][i]-W[k][j][im]/rho[k][j][im])+
                //idz*(U[k][j][i]/rho[k][j][i]-U[k-1][j][i]/rho[k-1][j][i])     
                     );
}
real_t sigma_zz(indx_t i,indx_t j,indx_t k)
{
    
    real_t ***rho =fields[0];
    real_t ***W   =fields[2];

    return 
    -pressure(rho[k][j][i])
    +gamma*lap(0,i,j,k)
    +2*eta*idz*(W[k][j][i]-W[k-1][j][i]);
   // +2*eta*idz*(W[k][j][i]/rho[k][j][i]-W[k-1][j][i]/rho[k-1][j][i]);
   
}
void computeDFieldsAtZ(indx_t k)
{
	indx_t i,j;
    real_t ***rho =fields[0];
    real_t ***U   =fields[1];
    real_t ***W   =fields[2];
    real_t ***pres=fields[3];
    
    real_t ***drho=dfields[0];
    real_t ***dU  =dfields[1];
    real_t ***dW  =dfields[2];
    real_t ***dpres=fields[3];
    
    
    
    indx_t km=k-1; /* always valid for there are ghosts */
    indx_t kp=k+1; /* always valid for there are ghosts */
 
    for(j=ymax;j>=ymin;--j)
    {
        indx_t jm = j-1;
        indx_t jp = j+1;
        if( jm < ymin ) jm = ymax;
        if( jp > ymax ) jp = ymin;
        for(i=xmax;i>=xmin;--i)
        {
            indx_t im = i-1;
            indx_t ip = i+1;
            if( im < xmin ) im = xmax;
            if( ip > xmax ) ip = xmin;
            
            drho[k][j][i]=-idx*(U[k][j][i]-U[k][j][im])-idz*(W[k][j][i]-W[km][j][i]);
            
            
            dU[k][j][i]  =
            idx*(sigma_xx(ip,j,k)-sigma_xx(i,j,k))+
            idz*0.5*(
                     (sigma_xz(i,j,kp)+sigma_xz(ip,j,kp))-
                     (sigma_xz(i,j,km)+sigma_xz(ip,j,km))
                     );
            
            dW[k][j][i]  =
            idz*(sigma_zz(i,j,kp)-sigma_zz(i,j,k))+
            idx*0.5*(
                     (sigma_xz(ip,j,kp)+sigma_xz(ip,j,k))-
                     (sigma_xz(im,j,kp)+sigma_xz(im,j,k))
                     );

 
            dpres[k][j][i]=0;
            pres[k][j][i]=pressure(rho[k][j][i]);
            
        }
        
    }
	
}

void computeDFieldsAtZOld(indx_t k)
{
	indx_t i,j;
    real_t ***rho =fields[0];
    real_t ***U   =fields[1];
    real_t ***W   =fields[2];
    real_t ***pres=fields[3];
 
    
    
    
    real_t ***drho=dfields[0];
    real_t ***dU  =dfields[1];
    real_t ***dW  =dfields[2];
    real_t ***dpres=fields[3];

    
 
    
    indx_t km=k-1; /* always valid for there are ghosts */
    indx_t kp=k+1; /* always valid for there are ghosts */
    for(j=ymax;j>=ymin;--j)
    {
        indx_t jm = j-1;
        indx_t jp = j+1;
        if( jm < ymin ) jm = ymax;
        if( jp > ymax ) jp = ymin;
        for(i=xmax;i>=xmin;--i)
        {
            indx_t im = i-1;
            indx_t ip = i+1;
            if( im < xmin ) im = xmax;
            if( ip > xmax ) ip = xmin;
            
            drho[k][j][i]=-idx*(U[k][j][i]-U[k][j][im])-idz*(W[k][j][i]-W[km][j][i]);
            
            
            dU[k][j][i]  =
            -idx*(pressure(rho[k][j][ip])-pressure(rho[k][j][i]))
            +idx*gamma*(lap(0,ip,j,k)-lap(0,i,j,k))
            +eta*lap(1,i,j,k)
            +noise*alea;
            dW[k][j][i]  =
            -idz*(pressure(rho[kp][j][i])-pressure(rho[k][j][i]))
            +idz*gamma*(lap(0,i,j,kp)-lap(0,i,j,k))
            +eta*lap(2,i,j,k)
            +noise*alea;
            

            dpres[k][j][i]=0;
            pres[k][j][i]=pressure(rho[k][j][i]);
            
        }
        
    }
	
}

void computeDFields(int bulk)
{
	indx_t k;
    
    if(bulk)
    {
        
        for(k=zmax-NG;k>=zmin+NG;--k)
            computeDFieldsAtZ(k);
    }
    else
    {
        for(k=0;k<NG;k++)
        {
            computeDFieldsAtZ(zmax-k);
            computeDFieldsAtZ(zmin+k);
        }
    }
    
}
static void integrate()
{
    
	size_t i;
	size_t j;
    
    for( i=0; i < NC; ++i )
        sendRequests(i);
    computeDFields(1); //bulk;

    for( i=0; i < NC; ++i )
        waitRequests(i);
    
    computeDFields(0); //boundaries;
    
    
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
}
double  mesureTimeForExchangingGhost()
{
	double elapsedTime;
	_BARRIER;
	elapsedTime=-MPI_Wtime();
	sendRequests(0);
	waitRequests(0);
	_BARRIER;
	elapsedTime+=MPI_Wtime();
	
	return elapsedTime;
    
}
static void initMetaBubble()
{
	indx_t i,j,k;
    real_t x,y,z;
    real_t ***rho=fields[0];
    real_t ***U  =fields[1];
    real_t ***V  =fields[2];
    srand(time(NULL)+rank);
	
    for(k=zmax;k>=zmin;--k)  
	{
        z=k*dz-Lz*0.5;
		for(j=ymax;j>=ymin;--j)
		{
            y=j*dy-Ly*0.5;
			for(i=xmax;i>=xmin;--i)
			{
                x=i*dx-Lx*0.5;
                if(x*x+z*z<(Lx*Lx+Lz*Lz)*0.01)
                    rho[k][j][i] =0.9+0.001*alea;
                else
                    rho[k][j][i] =1.1;

                U[k][j][i] =0;
                V[k][j][i] =0;
            }
		}
	}
}
#define initBubble 0
static void init_fields()
{
	indx_t i,j,k;
    real_t x,y,z;
    real_t ***rho=fields[0];
    real_t ***U  =fields[1];
    real_t ***V  =fields[2];
    srand(time(NULL)+rank);
	
    for(k=zmax;k>=zmin;--k)  
	{
        z=k*dz-Lz*0.5;
		for(j=ymax;j>=ymin;--j)
		{
            y=j*dy-Ly*0.5;
			for(i=xmax;i>=xmin;--i)
			{
                x=i*dx-Lx*0.5;
              //  if(x*x+z*z<100)
                if(initBubble)
                {
                    if(x*x+z*z<(Lx*Lx+Lz*Lz)*0.02)
                        rho[k][j][i] =0.2+0.00*alea;
                    else
                        rho[k][j][i] =1.1;
                    
                   // rho[k][j][i]=1.09+0.001*alea;
                }
                else
                {
                    rho[k][j][i] =0.6+0.000*alea;
                    U[k][j][i] =+0.001*alea;
                    V[k][j][i] =+0.001*alea;
                }
               
            }
		}
	}
}



void initSimulation(void)
{
    /***************************************************************************
	 * MPI setup
	 **************************************************************************/
	
	_CHECK(MPI_Comm_size(MPI_COMM_WORLD,&size) );
	_CHECK(MPI_Comm_rank(MPI_COMM_WORLD,&rank) );
	above = rank + 1; if( above >= size ) above = 0;
	below = rank - 1; if( below <  0    ) below = size-1;
	_BARRIER;
	fprintf( stderr, "-- ready %d.%d (%d->%d->%d)\n", rank, size, below, rank, above );
	fflush( stderr );
	
    
    num_reqs = 4* NC;
    
	requests = (MPI_Request *)calloc(num_reqs,sizeof(MPI_Request));
	/***************************************************************************
	 * Slicing along z, depending on rank and size
	 **************************************************************************/
	{
		indx_t left = Nz;         /* slices left  */
		indx_t from = 0;          /* first offset */
		indx_t todo = left/size;  /* first height */
		indx_t r;
		for(r=1;r<=rank;++r)
		{
			assert(r<size);
			from += todo;          /* forward offset */
			left -= todo;          /* decrease count */
			todo  = left/(size-r); /* next height    */
		}
		zmin = from;
		zmax = zmin + todo - 1;
		items_per_field = items_per_slice * todo;
		_BARRIER;
		fprintf(stderr,"-- for %2d/%2d: zmin=%4d, zmax=%4d\n", rank, size, (int)zmin, (int)zmax );
		fflush( stderr );
	}
	
	/***************************************************************************
	 * Allocating arrays with ghosts
	 **************************************************************************/
	zlo = zmin - NG;
	zhi = zmax + NG;
	if(  atexit(delete_fields) != 0 )
	{
		perror("atexit delete_fields");
		exit(-1);
	}
	create_fields();
	
	/***************************************************************************
	 * prepare requests
	 **************************************************************************/
	if( atexit(delete_requests) != 0 )
	{
		perror("atexit delete_requests");
		exit(-1);
	}
	create_requests();
	
	_BARRIER;
	if(rank==0) 
	{
		fprintf( stderr, "-- ready to compute\n");
		fflush(stderr);
	}
	
	
	/***************************************************************************
	 * initialize fields, initially at 0
	 **************************************************************************/
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
    for(i=0;i<100;i++)
    {
        integrate();
    }
    
    if(sim->savingFiles==1)
    {
        /* tentive for saving using python
        char  export_cmd[1024];
        const char *export_format =
        "d = ExportDBAttributes()\n"
        "d.db_type = \"Silo\"\n"
        "d.dirname=\"results\"\n"
        "d.filename = \"file%d\"\n"
        "d.variables=(\"rho\",\"u\")\n"    // Tuple of variables to export. "default" means the plotted variable.
        "ExportDatabase(d)\n";
        sprintf(export_cmd,export_format,sim->cycle);
        
       
        VisItExecuteCommand(export_cmd);
        
        sleep(0.2);
        _BARRIER;
        */
        
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
    
    //  sleep(1);
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
	
	dx   = Lx / ( Nx-1 );
	dy   = Ly / ( Ny-1 );
	dz   = Lz / ( Nz-1 );
	idx2 = 1/(dx*dx); 
	idy2 = 1/(dy*dy);
	idz2 = 1/(dz*dz);
	srand( time(NULL) );
	setvbuf( stderr, NULL, 0, _IONBF );
    
    
	
	_CHECK(MPI_Init(&argc,&argv));
	/***************************************************************************
	 * init of MPI, slicing of the domains, etc....
	 **************************************************************************/
    initSimulation(); 
    
    
    /***************************************************************************
     * Initialize environment variables. 
	 **************************************************************************/
    rmesh_dims[0]=Nx;
  //  rmesh_dims[1]=Ny;
    rmesh_dims[1]=zmax-zmin+1+NG;
    rmesh_dims[2]=1;
    
    
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
        
        if(VisItInitializeSocketAndDumpSimFile("bubble0",
                                               "Parallel C prototype simulation connects to VisIt",
                                               argv[0], NULL, "bubble0.ui", NULL))
      
        {
            fprintf(stderr,"VisItInitializeSocketAndDumpSimFile ok\n");
        }
        else
        {
            fprintf(stderr,"problem with VisItInitializeSocketAndDumpSimFile \n");
        }
    }
	init_fields();
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
