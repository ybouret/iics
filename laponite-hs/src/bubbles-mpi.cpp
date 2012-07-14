#include "bubbles.hpp"

void Bubbles:: dispatch_all(mpi &MPI)
{
    
    //==========================================================================
    // broadcast #bubbles
    //==========================================================================
    size_t num_bubbles = 0;
    if( MPI.IsMaster )
    {
        num_bubbles = count();
    }
    else
    {
        none();
    }
    MPI.BcastAs<size_t>(num_bubbles, 0, MPI_COMM_WORLD);
    
}