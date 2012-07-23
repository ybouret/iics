#include "yocto/utest/run.hpp"
#include "../segmenter.hpp"


YOCTO_UNIT_TEST_IMPL(seg)
{
    
    //==========================================================================
    //
    // prepare the grid
    // 
    //==========================================================================
    const Vertex  box(10,10);
    Bubbles       bubbles( box );
    const PBC    &pbc = bubbles.pbc;
    const Layout  lay( Coord(0,0), Coord(20,30) );
    const Region  reg( Vertex(0,pbc.lo), Vertex(box.x,pbc.up) );
    GhostsSetup   gs;
    FieldsSetup   fs;
    
    
}
YOCTO_UNIT_TEST_DONE()

