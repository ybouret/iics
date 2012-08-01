#include "bubble.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
void Bubble:: dispatch_topology( const mpi &MPI )
{
    if( MPI.IsParallel )
    {
        const bool master     = MPI.IsMaster;
        size_t     num_tracer = 0;
        
        //--------------------------------------------------------------------------
        // initialize
        //--------------------------------------------------------------------------
        if( master )
        {
            num_tracer = size;
        }
        else 
        {
            empty();
        }
        
        //--------------------------------------------------------------------------
        // broadcast #tracers
        //--------------------------------------------------------------------------
        MPI.__Bcast<size_t>(num_tracer, 0, MPI_COMM_WORLD);
        
        //--------------------------------------------------------------------------
        // broadcast physics
        //--------------------------------------------------------------------------
        MPI.Bcast(&pressure, IO_COUNT, MPI_REAL_TYPE, 0, MPI_COMM_WORLD);

        if(!master)
            append(num_tracer);
        
        //--------------------------------------------------------------------------
        // broadcast topology
        //--------------------------------------------------------------------------
        Tracer *p = root;
        for( size_t i=num_tracer;i>0;--i,p=p->next)
        {
            MPI.Bcast( &p->vertex, Tracer::IO_COUNT, MPI_REAL_TYPE, 0, MPI_COMM_WORLD);
        }
    }
    
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
void Bubble:: assemble_topology( const mpi &MPI )
{
    static const int tag = 1000;
    if( MPI.IsParallel )
    {
        const bool master   = MPI.IsMaster;
        if( master )
        {
            MPI_Status status;
            fprintf( stderr, "master changed=%u\n", unsigned(spots.size) );
            for( int r = 1; r < MPI.CommWorldSize; ++r )
            {
                //--------------------------------------------------------------
                // receive the num_changed
                //--------------------------------------------------------------
                
                const size_t num_changed = MPI.__Recv<size_t>(r, tag, MPI_COMM_WORLD, status);
                fprintf( stderr, "rank %d changed=%u\n", r, unsigned(num_changed) );
                
                //--------------------------------------------------------------
                // far all changed
                //--------------------------------------------------------------
                Tracer *p = root;
                for( size_t i=0; i < num_changed; ++i )
                {
                    //----------------------------------------------------------
                    // receive the #jump to perform
                    //----------------------------------------------------------
                    size_t jump = MPI.__Recv<size_t>(r, tag, MPI_COMM_WORLD, status);
                    
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
            const size_t num_changed = spots.size;
            MPI.__Send<size_t>(num_changed, 0, tag, MPI_COMM_WORLD);
            
            //------------------------------------------------------------------
            // far all changed
            //------------------------------------------------------------------
            const Spot *spot = spots.head;
            for( size_t i=0; i<num_changed; ++i , spot=spot->next)
            {
                //--------------------------------------------------------------
                // send the #jump to perform
                //--------------------------------------------------------------
                MPI.__Send<size_t>(spot->jump, 0, tag, MPI_COMM_WORLD);
                
                //--------------------------------------------------------------
                // send the new coordinates
                //--------------------------------------------------------------
                MPI.Send(& spot->handle->vertex, 2, MPI_REAL_TYPE, 0, tag, MPI_COMM_WORLD);
            }
        }
    }
    
}

