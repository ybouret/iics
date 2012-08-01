#include "yocto/utest/run.hpp"
#include "../cell.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/auto-ptr.hpp"
#include "yocto/swamp/mpi.hpp"

static inline double vproc( const double &x ) { return x; }

YOCTO_UNIT_TEST_IMPL(move)
{
    
    const mpi &MPI = mpi::init( &argc, &argv);
    
    Vertex box(10,10);
    Vertex center(box.x/2,0);
    
    Cell cell(20,30,box,MPI);
    if( cell.sim_master )
    {
        cell.bubbles.create()->map_peanut(center, 3.5, 0.9 + Alea() * 0.05 );
    }
    //SaveGrid( cell.mesh, vformat("grid%d.%d.dat", cell.sim_size,cell.sim_rank));
    
    const Layout B_layout( cell.sim_lower, Coord(cell.sim_upper.x,cell.sim_upper.y+1) );
    auto_ptr<Array> pB0;
    if( cell.sim_master )
        pB0.reset( new standalone<Array>(B_layout) );
    
    for( int iter=1; iter <= 100; ++iter )
    {
        cell.dispatch_all();
        if( cell.sim_master ) pB0->ldz();
        _mpi::collect0(MPI, pB0.__get(), cell.B, B_layout);
        if( cell.sim_master )
        {
            pB0->ppm( vformat("b%03d.ppm",iter), vformat("B"), B_layout, vproc);
        }
        
        
        Bubble *bubble = cell.bubbles.first();
        bubble->translate( Vertex(0,-0.1));
    }
    
    //cell.bubbles.first()->save_spots( vformat("spots%d.%d.dat",cell.sim_size,cell.sim_rank ) );
    //cell.segmenter.save_junctions( vformat("junc%d.%d.dat",cell.sim_size,cell.sim_rank) );
    //cell.bubbles.first()->save_inside( vformat("inside%d.%d.dat",cell.sim_size,cell.sim_rank), cell.mesh);
}
YOCTO_UNIT_TEST_DONE()
