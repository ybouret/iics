#include "../bubbles.hpp"

void Bubbles:: dispatch(const mpi &MPI )
{
    
    //--------------------------------------------------------------------------
    // broadcast num bubbles
    //--------------------------------------------------------------------------
    size_t num_bubbles = 0;
    if( MPI.IsFirst )
    {
        num_bubbles = bubbles.size;
    }
    MPI.Bcast(num_bubbles, 0, MPI_COMM_WORLD);
    
    //--------------------------------------------------------------------------
    // prepare bubbles on slaves
    //--------------------------------------------------------------------------
    if( !MPI.IsFirst )
    {
        create(num_bubbles);
    }
    
    assert( bubbles.size == num_bubbles );
    
    //--------------------------------------------------------------------------
    // dispatch each bubble
    //--------------------------------------------------------------------------
    for( Bubble *bubble = bubbles.head; bubble; bubble=bubble->next )
    {
        bubble->dispatch(MPI);
    }
    
    
}


void Bubbles:: assemble( const mpi & MPI )
{
#if !defined(NDEBUG)
    size_t num_bubbles = 0;
    if( MPI.IsFirst )
    {
        num_bubbles = bubbles.size;
    }
    MPI.Bcast(num_bubbles, 0, MPI_COMM_WORLD);
    assert( num_bubbles == bubbles.size );
#endif
    
    //--------------------------------------------------------------------------
    // assemble each bubble
    //--------------------------------------------------------------------------
    for( Bubble *bubble = bubbles.head; bubble; bubble=bubble->next )
    {
        bubble->assemble(MPI);
    }
    
    
}

