#include "yocto/utest/run.hpp"
#include "yocto/code/rand.hpp"
#include "shape.hpp"



static inline
void check( Bubble *b, const char *pfx )
{
    Shape::Rotate(b, numeric<Real>::pi * alea<Real>() );
    b->save_dat( vformat("%s0.dat", pfx ) );
    std::cerr << "Original Area=" << b->area << std::endl;
    b->adjust_contour();
    std::cerr << "Adjusted Area=" << b->area << std::endl;

    b->save_dat( vformat("%s1.dat", pfx ) );
    
    Shape::Grow(b, 2);
    std::cerr << "x2 Original Area=" << b->area << std::endl;
    b->save_dat( vformat("%s2.dat", pfx ) );
    b->adjust_contour();
    std::cerr << "x2 Adjusted Area=" << b->area << std::endl;
    b->save_dat( vformat("%s3.dat", pfx ) );

    Shape::Grow(b, 0.25);
    std::cerr << "/4 Original Area=" << b->area << std::endl;
    b->save_dat( vformat("%s4.dat", pfx ) );
    b->adjust_contour();
    std::cerr << "/4 Adjusted Area=" << b->area << std::endl;
    b->save_dat( vformat("%s5.dat", pfx ) );
    
    std::cerr << std::endl;
}

YOCTO_UNIT_TEST_IMPL(ctrl)
{
    Real   lam = 0.5;
    Real   gam = 1;
    Bubble bubble(lam,gam,1);
    
    std::cerr << "Circle" << std::endl;
    Shape::Circle(&bubble, Vertex(0,0), 4);
    check(&bubble, "c");
    
    std::cerr << "Ellipse" << std::endl;
    Shape::Ellipse(&bubble, Vertex(0,0), Vertex(5,2));
    check(&bubble,"e");
    
    std::cerr << "Blob" << std::endl;
    Shape::Blob(&bubble, Vertex(0,0), 4, 0.5 + 0.45 * alea<Real>(), 0.95 * alea<Real>() );
    check(&bubble,"b");

}
YOCTO_UNIT_TEST_DONE()
