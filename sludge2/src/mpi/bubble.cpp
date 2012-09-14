#include "../bubble.hpp"

void Bubble:: dispatch(const mpi &MPI)
{
    //--------------------------------------------------------------------------
    // empty spots
    //--------------------------------------------------------------------------
    spots.empty();
    
    //--------------------------------------------------------------------------
    // broadcast num tracers
    //--------------------------------------------------------------------------
    size_t num_tracers = 0;
    if( MPI.IsMaster )
    {
        num_tracers = size;
    }
    else
    {
        empty();
    }
    MPI.__Bcast(num_tracers, 0, MPI_COMM_WORLD );
    
    //--------------------------------------------------------------------------
    // create tracers on slave
    //--------------------------------------------------------------------------
    if( !MPI.IsMaster )
    {
        append( num_tracers );
    }
    
    //--------------------------------------------------------------------------
    // brodacast tracers
    //--------------------------------------------------------------------------
    Tracer *tracer = root;
    for( size_t i=num_tracers;i>0;--i,tracer=tracer->next)
    {
        MPI.Bcast( &tracer->vertex, Tracer::IO_COUNT, MPI_REAL_TYPE, 0, MPI_COMM_WORLD);
    }
    
}

void Bubble:: assemble( const mpi &MPI )
{
    static const int tag = 0xA55E;
#if !defined(NDEBUG)
    size_t num_tracers = 0;
    if( MPI.IsMaster )
        num_tracers = size;
    MPI.__Bcast(num_tracers, 0, MPI_COMM_WORLD );
    assert( size == num_tracers );
#endif
    
    if( MPI.IsMaster )
    {
        //----------------------------------------------------------------------
        // loop over slaves
        //----------------------------------------------------------------------
        MPI_Status status;
        for( int r = 1; r < MPI.CommWorldSize; ++r )
        {
            const size_t num_spots = MPI.__Recv<size_t>(r, tag, MPI_COMM_WORLD, status);
            fprintf( stderr, "from %d: #spots=%u\n", r, unsigned(num_spots));
            
            Tracer *tracer = root;
            for( size_t i=0; i < num_spots; ++i )
            {
                //--------------------------------------------------------------
                // get jump
                //--------------------------------------------------------------
                const size_t jump = MPI.__Recv<size_t>(r,tag,MPI_COMM_WORLD,status);
                
                //--------------------------------------------------------------
                // find right tracer
                //--------------------------------------------------------------
                for( size_t j=jump;j>0;--j) tracer = tracer->next;
                
                //--------------------------------------------------------------
                // recv vertex
                //--------------------------------------------------------------
                MPI.Recv( &(tracer->vertex), 2, MPI_REAL_TYPE, r, tag, MPI_COMM_WORLD, status);
            }
        }
    }
    else
    {
        //----------------------------------------------------------------------
        // send #slots
        //----------------------------------------------------------------------
        MPI.__Send<size_t>(spots.size, 0, tag, MPI_COMM_WORLD);
        for( const Spot *spot = spots.head; spot; spot=spot->next )
        {
            //------------------------------------------------------------------
            // send jump
            //------------------------------------------------------------------
            MPI.__Send<size_t>( spot->jump, 0, tag, MPI_COMM_WORLD);
            
            //------------------------------------------------------------------
            // send vertex
            //------------------------------------------------------------------
            MPI.Send( & (spot->handle->vertex), 2, MPI_REAL_TYPE, 0, tag, MPI_COMM_WORLD);
        }
        
    }
}