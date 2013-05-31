#include "yocto/utest/run.hpp"
#include "yocto/code/rand.hpp"
#include "junctions.hpp"
#include "shape.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/code/rand.hpp"
#include "yocto/spade/data-block.hpp"

YOCTO_UNIT_TEST_IMPL(segment)
{
    
    Real radius = 5.1;
    if( argc > 1 )
    {
        radius = strconv::to<Real>( argv[1], "radius" );
    }
    
    //==========================================================================
    //
    // Create the grid
    //
    //==========================================================================
    std::cerr << "*** Building grid" << std::endl;
    array_db       adb;
    const Layout   L( Coord(1,1), Coord(21,40) );
    Grid           grid(L,adb);
    const Region2D R( Vertex(-5,-5), Vertex(5,5) );
    
    Junctions junctions(grid);
    grid.regular_map_to(R,L);
    
    __Grid::SaveDat(grid, "grid.dat" );


    //==========================================================================
    //
    // Create the bubbles
    //
    //==========================================================================
    Bubbles bubbles;
    bubbles.lambda = __Grid::ComputeLambda(grid);
    std::cerr << "lambda=" << bubbles.lambda << std::endl;

    {
        const Vertex center(-2.5,-2.5);
        Bubble *bubble = bubbles.append();
        Shape::Blob(bubble, center, 3, 0.5+0.45*alea<Real>(), alea<Real>() );
        Shape::Rotate(bubble, numeric<Real>::two_pi * alea<Real>() );
    }
    
    {
        const Vertex center(2.5,2.5);
        Bubble *bubble = bubbles.append();
        Shape::Blob(bubble, center, 3, 0.5+0.45*alea<Real>(), alea<Real>() );
        Shape::Rotate(bubble, numeric<Real>::two_pi * alea<Real>() );
    }
    
    bubbles.regularize_all();
    bubbles.collect_all_markers(grid.Y()[grid.lower.y], grid.Y()[grid.upper.y]);
    
    for(const Bubble *b = bubbles.head; b; b=b->next)
    {
        std::cerr << "Bubble #" << b->UID << std::endl;
        const string pfx = vformat("b%u", unsigned(b->UID));
        b->save_all(pfx);
    }
    
    
    //==========================================================================
    //
    // Segment
    //
    //==========================================================================
    junctions.load(bubbles);
    junctions.save_dat("j.dat");
    junctions.save_t("j_t.dat");
    junctions.save_n("j_n.dat");
    
    standalone<Array> B(L);
    
    junctions.segment(B);
    junctions.save_inside_of(B, "inside.dat");
    
    std::cerr << "sizeof(Tracer)   = " << sizeof(Tracer)   << std::endl;
    std::cerr << "sizeof(Junction) = " << sizeof(Junction) << std::endl;
    std::cerr << "sizeof(Marker)   = " << sizeof(Marker)   << std::endl;
    std::cerr << "sizeof(Bubble)   = " << sizeof(Bubble)   << std::endl;
    
}
YOCTO_UNIT_TEST_DONE()
