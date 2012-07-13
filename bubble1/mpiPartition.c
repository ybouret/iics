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

void initSimulation(void)
{ 
    int i;
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
    
    
    rmesh_dims[0]=Nx;
    rmesh_dims[1]=zmax-zmin+1+NG;
    rmesh_dims[2]=1;
    for(i=0;i<3;i++)
        rmesh[i] = (float *)malloc(sizeof(float) * rmesh_dims[i]);
    
    for(i = 0; i < rmesh_dims[0]; ++i)
        rmesh[0][i] = i*dx;
    for(i = 0; i < rmesh_dims[1]; ++i)
        rmesh[1][i] = (zmin+i)*dz;
}