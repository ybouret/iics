#include "yocto/utest/run.hpp"
#include "../bubble.hpp"
#include "yocto/ios/ocstream.hpp"

static void save_spots( const Bubble &bubble, const string &filename )
{
    ios::ocstream fp(filename,false);
    for( const Spot *spot = bubble.spots.head; spot; spot=spot->next )
    {
        fp("%g %g\n", spot->handle->vertex.x, spot->handle->vertex.y);
    }
    
}

YOCTO_UNIT_TEST_IMPL(bubble)
{
    
    const mpi &MPI = mpi::init( &argc, &argv);
    
    const int rank = MPI.CommWorldRank;
    const int size = MPI.CommWorldSize;
    
    double lambda = 1;
    Vertex box(10,10);
    PBC    pbc(box.y);
    Tracer::Cache tcache;
    Spot::Cache   scache;
    
    Vertex center( box.x/2, 0.0 );
    
    Bubble bubble(lambda,pbc,tcache,scache);
    
    if( MPI.IsMaster )
    {
        bubble.map_circle( center, box.x/3);
    }
    
    
    const Real y_lo = pbc.lo + (rank*pbc.L)/size;
    const Real y_up = pbc.lo + ((rank+1)*pbc.L)/size;
    MPI.Printf(stderr, "rank %d> %g -> %g\n",rank, y_lo, y_up);
    
    bubble.dispatch_topology(MPI);
    bubble.mark_and_find_spots_within(y_lo, y_up);
    bubble.compute_geometry();
    MPI.Printf(stderr,"rank %d> #spots= %lu\n", rank, bubble.spots.size );
    
    bubble.save_dat( vformat("bubble%d.%d.dat",rank,size) );
    save_spots(bubble, vformat("spots%d.%d.dat",rank,size) );
    
    bubble.assemble_topology(MPI);
    
}
YOCTO_UNIT_TEST_DONE()

