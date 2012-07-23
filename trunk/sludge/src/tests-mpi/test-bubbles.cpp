#include "yocto/utest/run.hpp"
#include "../bubbles.hpp"
#include "yocto/ios/ocstream.hpp"

YOCTO_UNIT_TEST_IMPL(bubbles)
{
    
    const mpi &MPI = mpi::init( &argc, &argv);
    
    const int rank = MPI.CommWorldRank;
    const int size = MPI.CommWorldSize;
    
    Vertex box(10,10);
    Vertex center(box.x/2,0);
    Bubbles bubbles(box);
    const PBC &pbc = bubbles.pbc;
    if( MPI.IsMaster )
    {
        bubbles.create()->map_circle( center, 1.4 );
        bubbles.create()->map_circle( center, 3   );
    }

    bubbles.check_and_dispatch_all(MPI);
    
    const Real y_lo = pbc.lo + (rank*pbc.L)/size;
    const Real y_up = pbc.lo + ((rank+1)*pbc.L)/size;
    MPI.Printf(stderr, "rank %d> %g -> %g\n",rank, y_lo, y_up);
    bubbles.check_geometries_within(y_lo, y_up);
    
    for( Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next )
    {
        for( Spot *spot = bubble->spots.head; spot; spot=spot->next )
        {
            Vertex &v = spot->handle->vertex;
            v.x += 0.1 * (Alea() - 0.5 );
            v.y += 0.1 * (Alea() - 0.5 );
        }
        
    }
    bubbles.assemble_all(MPI);
    
    
}
YOCTO_UNIT_TEST_DONE()
