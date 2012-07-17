#include "bubbles.hpp"

void Bubbles:: dispatch_all(const mpi &MPI)
{
    const bool master   = MPI.IsMaster;
    //==========================================================================
    // broadcast #bubbles
    //==========================================================================
    size_t num_bubbles = 0;
    if( master )
    {
        num_bubbles = count();
    }
    MPI.BcastAs<size_t>(num_bubbles, 0, MPI_COMM_WORLD);
    
    //==========================================================================
    // adjust resources
    //==========================================================================
    if( !master )
    {
        none();
        create(num_bubbles);
    }
    
    //==========================================================================
    // broadcats all bubbles
    //==========================================================================
    for( Bubble *b = b_list.head;b;b=b->next)
    {
        b->dispatch(MPI);
    }
    
    
}

void Bubbles:: assemble_all(const mpi &MPI )
{
#if !defined (NDEBUG)
    const bool master   = MPI.IsMaster;
    //==========================================================================
    // broadcast #bubbles
    //==========================================================================
    size_t num_bubbles = 0;
    if( master )
    {
        num_bubbles = count();
    }
    MPI.BcastAs<size_t>(num_bubbles, 0, MPI_COMM_WORLD);
    if( num_bubbles != count() )
        throw exception("Invalid bubbles count to collect rank %d", MPI.CommWorldRank);
#endif
    
    for( Bubble *b = b_list.head;b;b=b->next)
    {
        b->assemble(MPI);
    }
    
}
