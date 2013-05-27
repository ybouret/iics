#include "yocto/utest/run.hpp"
#include "mpi/bubbles.hpp"

YOCTO_UNIT_TEST_IMPL(bubbles)
{
    YOCTO_MPI;
    MPI.Printf(stderr,"Starting...\n");
    Real lambda = 0.5;
    Bubbles bubbles(lambda);
    if( MPI.IsFirst )
    {
        
    }
    
    Parallel::BubblesBcast(MPI, bubbles);
    Hasher h;
    h.set();
    bubbles.hash(h);
    MPI.Printf(stderr, "Key=%d\n", h.getKey());
}
YOCTO_UNIT_TEST_DONE()

