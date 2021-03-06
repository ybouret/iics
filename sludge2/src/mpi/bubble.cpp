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
    if( MPI.IsFirst )
    {
        num_tracers = size; // from master
    }
    else
    {
        empty();            // for slaves
    }
    MPI.Bcast(num_tracers, 0, MPI_COMM_WORLD );
    
    //--------------------------------------------------------------------------
    // create tracers on slave
    //--------------------------------------------------------------------------
    if( !MPI.IsFirst )
    {
        for(size_t i=num_tracers;i>0;--i)
            append()->bubble = this;
    }
    
    //--------------------------------------------------------------------------
    // brodacast tracers
    //--------------------------------------------------------------------------
    Tracer *tracer = root;
    for( size_t i=num_tracers;i>0;--i,tracer=tracer->next)
    {
        assert(tracer->bubble==this);
        MPI.Bcast( &tracer->vertex, Tracer::IO_COUNT, MPI_REAL_TYPE, 0, MPI_COMM_WORLD);
    }
    
    //--------------------------------------------------------------------------
    // brodacast pressure sqq..
    //--------------------------------------------------------------------------
    MPI.Bcast(&pressure, Bubble::IO_COUNT, MPI_REAL_TYPE, 0, MPI_COMM_WORLD);
    
}

void Bubble:: assemble( const mpi &MPI )
{
    static const int tag = 0xA55E;
#if !defined(NDEBUG)
    size_t num_tracers = 0;
    if( MPI.IsFirst )
        num_tracers = size;
    MPI.Bcast(num_tracers, 0, MPI_COMM_WORLD );
    assert( size == num_tracers );
#endif
    
    MPI.Printf0( stderr, "\t[BEFORE] #%u: area=%g, pressure=%g\n", id, area, pressure);

    if( MPI.IsFirst )
    {
        //----------------------------------------------------------------------
        // loop over slaves
        //----------------------------------------------------------------------
        MPI_Status status;
        for( int r = 1; r < MPI.CommWorldSize; ++r )
        {
            const size_t num_spots = MPI.Recv<size_t>(r, tag, MPI_COMM_WORLD, status);
            Tracer *tracer = root;
            for( size_t i=0; i < num_spots; ++i )
            {
                //--------------------------------------------------------------
                // get jump
                //--------------------------------------------------------------
                const size_t jump = MPI.Recv<size_t>(r,tag,MPI_COMM_WORLD,status);
                
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
        //----------------------------------------------------------------------
        // TODO: update content
        //----------------------------------------------------------------------
        
        
        //----------------------------------------------------------------------
        // update area, iso content
        //----------------------------------------------------------------------
        update_area_full();
        pressure = content / area;
    }
    else
    {
        //----------------------------------------------------------------------
        // send #slots
        //----------------------------------------------------------------------
        MPI.Send<size_t>(spots.size, 0, tag, MPI_COMM_WORLD);
        for( const Spot *spot = spots.head; spot; spot=spot->next )
        {
            //------------------------------------------------------------------
            // send jump
            //------------------------------------------------------------------
            MPI.Send<size_t>( spot->jump, 0, tag, MPI_COMM_WORLD);
            
            //------------------------------------------------------------------
            // send vertex
            //------------------------------------------------------------------
            MPI.Send( & (spot->handle->vertex), 2, MPI_REAL_TYPE, 0, tag, MPI_COMM_WORLD);
        }
        
    }
    MPI.Printf0( stderr, "\t[AFTER ] #%u: area=%g, pressure=%g\n", id, area, pressure);

}