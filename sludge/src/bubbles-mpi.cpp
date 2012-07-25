#include "bubbles.hpp"


void Bubbles:: check_and_dispatch_all(const mpi &MPI)
{
    const bool master = MPI.IsMaster;
    
    //==========================================================================
    // determine #bubbles
    //==========================================================================
    size_t num_bubbles = 0;
    if( master )
    {
        check_topologies();
        num_bubbles = count();
    }
    else 
    {
        empty();
    }
    
    MPI.__Bcast<size_t>(num_bubbles, 0, MPI_COMM_WORLD);
    
    //==========================================================================
    // dispatch each bubble
    //==========================================================================
    if(!master)
    {
        for(size_t i=num_bubbles;i>0;--i) (void) create();
    }
    
    for( Bubble *bubble = first(); bubble; bubble = bubble->next )
    {
        bubble->dispatch_topology(MPI);
    }
        
}


void Bubbles:: assemble_all( const mpi &MPI )
{
    
#if !defined (NDEBUG)
    size_t num_bubbles = 0;
    if(MPI.IsMaster)
        num_bubbles = count();
    MPI.__Bcast<size_t>(num_bubbles,0,MPI_COMM_WORLD);
    assert( count() == num_bubbles );
#endif
    
    for( Bubble *bubble = first(); bubble; bubble = bubble->next )
    {
        bubble->assemble_topology(MPI);
    }
    
}
