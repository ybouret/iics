#include "yocto/utest/run.hpp"
#include "../cell.hpp"

#include "yocto/ios/ocstream.hpp"

YOCTO_UNIT_TEST_IMPL(cell)
{
    
    const mpi &MPI = mpi::init( &argc, &argv);
        
    Vertex box(10,10);
    Vertex center(box.x/2,0);
    
    Cell cell(20,30,box,MPI);
    if( cell.sim_master )
    {
        cell.bubbles.create()->map_peanut(center+Vertex(0,-4), 3.5, 0.95);
    }
    
    cell.dispatch_all();
    
    SaveGrid( cell.mesh, vformat("grid%d.%d.dat", cell.sim_size,cell.sim_rank));
    
    cell.bubbles.first()->save_spots( vformat("spots%d.%d.dat",cell.sim_size,cell.sim_rank ) );
    cell.segmenter.save_junctions( vformat("junc%d.%d.dat",cell.sim_size,cell.sim_rank) );
    cell.bubbles.first()->save_inside( vformat("inside%d.%d.dat",cell.sim_size,cell.sim_rank), cell.mesh);
}
YOCTO_UNIT_TEST_DONE()
