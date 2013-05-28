#include "yocto/utest/run.hpp"
#include "yocto/code/rand.hpp"
#include "grid.hpp"


YOCTO_UNIT_TEST_IMPL(grid)
{
    array_db       adb;
    const Layout   L( Coord(1,1), Coord(10,20) );
    Grid           grid(L,adb);
    const Region2D R( Vertex(-5,5), Vertex(-6,6) );
    
    grid.regular_map_to(R,L);
    
    __Grid::SaveDat(grid, "grid.dat" );
    
}
YOCTO_UNIT_TEST_DONE()

