#include "yocto/utest/run.hpp"
#include "../segmenter.hpp"
#include "yocto/code/utils.hpp"

YOCTO_UNIT_TEST_IMPL(seg)
{
    const unit_t NX  = 10;
    unit_t       NY  = 20;
    while( NY & 1 ) ++NY;
    const unit_t HNY = NY>>1;
    const Real   LX  = 4;
    const Real   LY  = 6;
    const Real   dX  = LX/NX;
    const Real   dY  = LY/NY;
    const PBC    pbc(LY);
    const Layout L( Coord(0,-HNY), Coord(NX,HNY) );
    array_db     adb;
    Grid         G(L,adb);
    Segmenter    seg(G);
    
    for( unit_t i=G.lower.x;i<=G.upper.x;++i)
    {
        G.X()[i] = (i*LX)/NX;
    }
    
    for( unit_t j=G.lower.y;j<=G.upper.y;++j)
    {
        G.Y()[j] = (j*LY)/NY;
    }
    std::cerr << "seg.Y[lower]=" << seg.Y[seg.Y.lower] << " / " << pbc.lo << std::endl;
    std::cerr << "seg.Y[upper]=" << seg.Y[seg.Y.upper] << " / " << pbc.up << std::endl;

    const Vertex center( LX/2, 0 );
    const Real   radius=  min_of<Real>(LX,LY)/4;
    Bubbles      bubbles(pbc);
    bubbles.lambda = min_of<Real>(dX,dY)/2;
    Bubble      *bubble = bubbles.append();
    bubble->map_peanut(center, radius, 0.85);
    bubble->compute_contour();
    
    SaveGrid(G, "grid.dat" );
    bubble->save_dat( "b.dat");
    
    
}
YOCTO_UNIT_TEST_DONE()

