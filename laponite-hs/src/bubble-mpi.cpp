#include "bubble.hpp"


void Bubble:: dispatch(mpi &MPI)
{
    
    if( MPI.CommWorldSize >0 )
    {
        //======================================================================
        // broadcast the number of points
        //======================================================================
        const bool master   = MPI.IsMaster;
        unsigned num_points = 0;
        if( master )
        {
            num_points = size;
        }
        else 
        {
            empty();
        }
        MPI.Bcast(&num_points, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI.Printf( stderr, "Rank %d> num_points=%u\n", MPI.CommWorldRank, num_points );
        
        //======================================================================
        // broadcast the components
        //======================================================================
        Point *p = root;
        for( unsigned i=0; i < num_points; ++i )
        {
            if( !master )
            {
                p = create();
                push_back(p);
            }
            assert(p!=NULL);
            V2D &v = *p;
            
            MPI.Bcast(&v, 2, MPI_REAL_TYPE, 0, MPI_COMM_WORLD);
            
            if(master)
                p = p->next;
        }
        MPI.Printf(stderr,"Rank %d> bubble.size=%u\n", MPI.CommWorldRank, unsigned(size) );
    }
}


void Bubble:: collect(mpi &MPI)
{
    static const int tag = 1000;
    if( MPI.CommWorldSize >0 )
    {
        const bool master   = MPI.IsMaster;
        MPI.Barrier(MPI_COMM_WORLD);
        if( master )
        {
            MPI_Status status;
            fprintf( stderr, "master changed=%u\n", unsigned(spots.size) );
            for( int r = 1; r < MPI.CommWorldSize; ++r )
            {
                unsigned num_changed = 0;
                MPI.Recv(&num_changed, 1, MPI_UNSIGNED, r, tag, MPI_COMM_WORLD, status);
                fprintf( stderr, "rank %d changed=%u\n", r, num_changed );
            }
        }
        else 
        {
            const unsigned num_changed = spots.size;
            MPI.Send(&num_changed, 1, MPI_UNSIGNED, 0, tag, MPI_COMM_WORLD);
        }
    }
}
