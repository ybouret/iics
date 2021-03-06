#include <mpi.h>
#include "arrays.h"
/*** For Silo ***/
#include <silo.h>
#include <stdio.h>
/* Instrumenting with visit */






static const size_t NC = 1;     /*!< two components */
const char                *cpntName[] = {"u","v"};  /*Name of the components*/
static indx_t       Nx = 256;   /*!< 0 -> Nx-1      */
static indx_t       Ny = 256;   /*!< 0 -> Ny-1      */
static indx_t       Nz = 256;   /*!< 0 -> Nz-1      */
static indx_t       NG = 1;     /*!< #ghosts        */
static real_t       Lx = Nx*1.0;
static real_t       Ly = Ny*1.0;
static real_t       Lz = Nz*1.0;

static real_t   ****fields          = NULL;
static real_t    ***laplacian       = NULL;
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
static real_t       idx2,idy2,idz2;
static real_t       dt = 0.0;

#define _CHECK(MPI_ROUTINE) do { const int __rc = MPI_ROUTINE; if( MPI_SUCCESS != __rc ) { fprintf( stderr, "Failure: " #MPI_ROUTINE "\n" ); exit(-1); } } while(0)
#define _BARRIER _CHECK(MPI_Barrier(MPI_COMM_WORLD))
static MPI_Request *requests = NULL;
static size_t       num_reqs = 0;
static const int    diff_tag = 7;
real_t paramMu=-1;
real_t paramA=0;

int   rmesh_dims[3];


#include "stubs.c"
#include "ui.cpp"
#include "writeToSilo.c"




static void create_fields()
{
	/***************************************************************************
	 ** NC component + 1 field for the laplacian value
	 **************************************************************************/
	size_t i;
	fields = (real_t ****)calloc(NC+1,sizeof(real_t ***));
	if( !fields )
	{
		perror("fields level-1");
		exit(-1);
	}
	for(i=0;i<=NC;++i)
	{
		fields[i] = icp_create_array3D(xmin,xmax,ymin,ymax,zlo,zhi);
		if( !fields[i] )
		{
			while( i > 0 )
			{
				icp_delete_array3D(fields[--i],xmin,xmax,ymin,ymax,zlo,zhi);
			}
			free(fields);
			fields=NULL;
			perror("fields level-2");
			exit(-1);
		}
	}
	
	/** the laplacian is the extraneous field **/
	laplacian = fields[NC];
}

static void delete_fields()
{
	if( fields )
	{
		size_t i=NC+1;
		while( i > 0 )
		{
			icp_delete_array3D(fields[--i],xmin,xmax,ymin,ymax,zlo,zhi);
		}
		free(fields);
		fields = NULL;
	}
	laplacian = NULL;
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

/*
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
 */
/*****************************************************************************************************
 *     Here we start the send/recv requests for variable i
 *****************************************************************************************************/
static void sendRequests(int i)
{
	const size_t nitems = items_per_slice * NG;
    const size_t j = i * 4;
    
	if( !requests )
	{
		perror("requests in sendRequests");
		exit(-1);
	}
    
    // send information to below
    MPI_Isend( &fields[i][zmin][ymin][xmin],    nitems, ICP_REAL, below, diff_tag, MPI_COMM_WORLD, &requests[j+0] );
    
    // send information to above
    MPI_Isend( &fields[i][zmax+1-NG][ymin][xmin], nitems, ICP_REAL, above, diff_tag, MPI_COMM_WORLD, &requests[j+1] );
    
    // recv information from below
    MPI_Irecv( &fields[i][zlo][ymin][xmin],     nitems, ICP_REAL, below, diff_tag, MPI_COMM_WORLD, &requests[j+2] );
    
    // recv information from above
    MPI_Irecv( &fields[i][zmax+1][ymin][xmin],  nitems, ICP_REAL, above, diff_tag, MPI_COMM_WORLD, &requests[j+3] );
	
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



static void compute_laplacianAtZ( real_t ***f, indx_t k)
{
    indx_t i,j;
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
            
            const real_t f0  = f[k][j][i];
            const real_t tf0 = f0 + f0;
            const real_t lx    = idx2*(f[k][j][im]-tf0+f[k][j][ip]);
            const real_t ly    = idy2*(f[k][jm][i]-tf0+f[k][jp][i]);
            const real_t lz    = idz2*(f[km][j][i]-tf0+f[kp][j][i]);
            laplacian[k][j][i] = lx + ly + lz;
        }
    }
}
/**************************************************************************************************
 * we compute the laplacian of the field f. bulk=1 means we compute the bulk,
 **************************************************************************************************/
static void compute_laplacian2( real_t ***f, int bulk)
{
	indx_t k;
    
    if(bulk)
    {
        for(k=zmax-NG;k>=zmin+NG;--k)
            compute_laplacianAtZ(f,k);
    }
    else
    {
        for(k=0;k<NG;k++)
        {
            compute_laplacianAtZ(f,zmax-k);
            compute_laplacianAtZ(f,zmin+k);
        }
    }
    
}




//static void reaction(){}
/*
 static void diffusion()
 {
 size_t i;
 size_t j;
 
 exchange_ghosts();
 for( i=0; i < NC; ++i )
 {
 real_t ***f = fields[i];
 compute_laplacian(f);
 {
 real_t       *dst = &f[zmin][ymin][xmin];
 const real_t *src = &laplacian[zmin][ymin][xmin];
 for( j=0; j < items_per_field; ++j )
 {
 
 dst[j] += dt * (src[j]+dst[j]-dst[j]*dst[j]*dst[j]);
 }
 }
 }
 }
 */
static void diffusion2()
{
	size_t i;
	size_t j;
    
    
    for( i=0; i < NC; ++i )
	{
        real_t ***f = fields[i];
        sendRequests(i);
		compute_laplacian2(f,1); //bulk
        waitRequests(i);
        compute_laplacian2(f,0); //boundaries
        
        
        real_t       *dst = &f[zmin][ymin][xmin];
        const real_t *src = &laplacian[zmin][ymin][xmin];
        for( j=0; j < items_per_field; ++j )
        {
            dst[j] += dt * (src[j]+dst[j]*paramMu-dst[j]*dst[j]*dst[j]+paramA);
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

#define alea (rand()/(0.0+RAND_MAX)*2-1)

static void init_fields()
{
	indx_t i,j,k;
    real_t x,y,z;
    
    if(!rank)
    {
        fprintf( stderr, "\tinit_fields\n");
    }
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
                if(x*x+y*y+z*z<100)
                    fields[0][k][j][i] =1;
                else
                    fields[0][k][j][i] =alea;
            }
		}
	}
    sendRequests(0);
    waitRequests(0);
}



void initSimulation(void)
{
    /***************************************************************************
	 * MPI setup
	 **************************************************************************/
	
	_CHECK(MPI_Comm_size(MPI_COMM_WORLD,&size) );
	_CHECK(MPI_Comm_rank(MPI_COMM_WORLD,&rank) );
	above = rank + 1;
    if( above >= size )
        above = 0;
	below = rank - 1;
    if( below <  0    )
        below = size-1;
	_BARRIER;
	fprintf( stderr, "-- ready %d.%-2d (%2d->%2d->%2d)\n", size, rank, below, rank, above );
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
    double elapsedTime=0;
    ++sim->cycle;
    sim->time += dt;
    if(rank==0)
        elapsedTime=MPI_Wtime();
    
    // Diffusion
    
    for(i=0;i<1;i++)
    {
        diffusion2();
    }
    
    
    if(sim->savingFiles==1)
    {
        writeDomain(sim->cycle);
        write_master(sim->cycle);
    }
    
    
    if(sim->visitIsConnected==1)
    {
        VisItTimeStepChanged();
        VisItUpdatePlots();
    }
    
    if(rank==0)
    {
        const double ell = MPI_Wtime() - elapsedTime;
        const double spd = 1.0 / ell;
        fprintf(stderr,"%2d cores:simulating: cycle=%6d, time=%10.5lf in %10.5lf s | speed= %10.5f /s\n",sim->par_size, sim->cycle, sim->time, ell, spd);
        fflush(stderr);
    }
    
    //  sleep(1);
}

#define MIN_SIZE 2

int main(int argc, char *argv[] )
{
	//unsigned count = 0;
    //int toWrite=0;
    //double startTime,endTime;
    char *env = NULL;
    
    
    switch( argc )
    {
            
        case 2: //-- same Nx, Ny, Nz
        {
            const indx_t user_N = atoi(argv[1]);
            if(user_N<=MIN_SIZE)
            {
                fprintf( stderr, "Invalid Nx=%d\n", int(user_N) );
                return 1;
            }
            Nx = Ny = Nz = user_N;
        }
            break;
            
        case 4: //-- Nx, Ny, Nz
        {
            Nx = atoi( argv[1] );
            Ny = atoi( argv[2] );
            Nz = atoi( argv[3] );
            if( Nx <= MIN_SIZE || Ny <= MIN_SIZE || Nz <= MIN_SIZE )
            {
                fprintf( stderr, "Invalid Nx=%d,Ny=%d,Nz=%d\n", int(Nx), int(Ny), int(Nz) );
                return 1;
            }
        }
            
            //-- no args
        case 1:
            break;
    }
    
    Lx = double(Nx);
    Ly = double(Ny);
    Lz = double(Nz);
    
    
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
	dt=0.1;
	srand( time(NULL) );
	setvbuf( stderr, NULL, 0, _IONBF );
    
    
    
	_CHECK(MPI_Init(&argc,&argv));
    /*
     fprintf( stderr, "\n");
     fprintf( stderr, "\t----------------\n");
     fprintf( stderr, "\tStarting Program\n");
     fprintf( stderr, "\t----------------\n");
     fprintf( stderr, "\n");
     */
    fprintf( stderr, "-- Initializing with Nx=%d, Ny=%d, Nz=%d\n", int(Nx), int(Ny), int(Nz));
	/***************************************************************************
	 * init of MPI, slicing of the domains, etc....
	 **************************************************************************/
    initSimulation();
    
    
    /***************************************************************************
     * Initialize environment variables.
	 **************************************************************************/
    rmesh_dims[0]=Nx;
    rmesh_dims[1]=Ny;
    rmesh_dims[2]=zmax-zmin+1+NG;
    
    simulation_data sim;
    simulation_data_ctor(&sim,rank,size);
    //VisItSetupEnvironment();
    
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
    
    if(sim.par_rank == 0)
        env = VisItGetEnvironment();
    VisItSetupEnvironment2(env);
    if(env != NULL)
        free(env);
    
    if(rank==0)
    {
        
        if(VisItInitializeSocketAndDumpSimFile("v0",
                                               "Parallel C prototype simulation connects to VisIt",
                                               argv[0], NULL, "v0.ui", NULL))
            
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
    //    set_interface(&sim);
    
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
