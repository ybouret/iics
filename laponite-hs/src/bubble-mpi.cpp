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
            
            
            MPI.Bcast(& p->vertex, 2, MPI_REAL_TYPE, 0, MPI_COMM_WORLD);
            
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
                //--------------------------------------------------------------
                // receive the num_changed
                //--------------------------------------------------------------
                unsigned num_changed = 0;
                MPI.Recv(&num_changed, 1, MPI_UNSIGNED, r, tag, MPI_COMM_WORLD, status);
                fprintf( stderr, "rank %d changed=%u\n", r, num_changed );
                
                //--------------------------------------------------------------
                // far all changed
                //--------------------------------------------------------------
                Point *p = root;
                for( unsigned i=0; i < num_changed; ++i )
                {
                    //----------------------------------------------------------
                    // receive the #jump to perform
                    //----------------------------------------------------------
                    unsigned jump = 0;
                    MPI.Recv(&jump, 1, MPI_UNSIGNED, r, tag, MPI_COMM_WORLD, status);
                    
                    //----------------------------------------------------------
                    // move forward
                    //----------------------------------------------------------
                    while(jump-->0) 
                    { 
                        assert(p!=NULL);
                        p=p->next; 
                    }
                    assert(p!=NULL);
                    
                    //----------------------------------------------------------
                    // receive the new coordinates
                    //----------------------------------------------------------
                    MPI.Recv(& p->vertex, 2, MPI_REAL_TYPE, r, tag, MPI_COMM_WORLD, status);
                }
            }
        }
        else 
        {
            //------------------------------------------------------------------
            // send the num_changed
            //------------------------------------------------------------------
            const unsigned num_changed = spots.size;
            MPI.Send(&num_changed, 1, MPI_UNSIGNED, 0, tag, MPI_COMM_WORLD);
            
            //------------------------------------------------------------------
            // far all changed
            //------------------------------------------------------------------
            const Spot *spot = spots.head;
            for( unsigned i=0; i<num_changed; ++i , spot=spot->next)
            {
                //--------------------------------------------------------------
                // send the #jump to perform
                //--------------------------------------------------------------
                MPI.Send( &spot->jump, 1, MPI_UNSIGNED, 0, tag, MPI_COMM_WORLD);
                
                //--------------------------------------------------------------
                // send the new coordinates
                //--------------------------------------------------------------
                MPI.Send(& spot->point->vertex, 2, MPI_REAL_TYPE, 0, tag, MPI_COMM_WORLD);
            }
        }
    }
}
