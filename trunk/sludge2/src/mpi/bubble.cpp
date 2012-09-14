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