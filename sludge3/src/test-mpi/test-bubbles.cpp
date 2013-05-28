#include "yocto/utest/run.hpp"
#include "mpi/bubbles.hpp"
#include "shape.hpp"

YOCTO_UNIT_TEST_IMPL(bubbles)
{
    YOCTO_MPI;
    MPI.Printf(stderr,"Starting...\n");
    Real lambda = 0.5;
    
    // create the bubble(s)
    Bubbles bubbles(lambda);
    if( MPI.IsFirst )
    {
        Shape::Circle(bubbles.append(),Vertex(0,0),1.3);
        Shape::Circle(bubbles.append(),Vertex(2,0),1.7);
    }
    
    // broadcast
    ParallelBubbles:: Bcast(MPI, bubbles);
    
    for( const Bubble *b = bubbles.head; b; b=b->next )
    {
		b->save_dat( vformat("bubble%d.%d.dat", MPI.CommWorldSize, MPI.CommWorldRank) );
	}
    
    // check
    Hasher h;
    h.set();
    bubbles.hash(h);
    MPI.Printf(stderr, "Key=%d\n", h.getKey());
}
YOCTO_UNIT_TEST_DONE()

