#include "yocto/utest/run.hpp"
#include "mpi/bubbles.hpp"
#include "shape.hpp"

YOCTO_UNIT_TEST_IMPL(bubbles)
{
    YOCTO_MPI;
    MPI.Printf(stderr,"Starting...\n");
    
    // create the bubble(s)
    Bubbles bubbles;
    bubbles.lambda = 0.5;
    bubbles.gamma  = 0.1;
    
    if( MPI.IsFirst )
    {
        Shape::Circle(bubbles.append(),Vertex(0,0),1.3);
        Shape::Circle(bubbles.append(),Vertex(2,0),1.7);
        for( Bubble *b = bubbles.head; b; b=b->next )
        {
            b->regularize();
        }
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
    
    
    // test collect, meaningless in that case
    ParallelBubbles::Collect(MPI, bubbles);
    h.set();
    bubbles.hash(h);
    MPI.Printf0(stderr, "Key0=%d\n", h.getKey());
    
}
YOCTO_UNIT_TEST_DONE()

