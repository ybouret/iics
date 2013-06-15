#include "yocto/utest/run.hpp"
#include "yocto/code/rand.hpp"
#include "shape.hpp"



static inline
void check( Bubble *b, const char *pfx )
{
    b->save_dat( vformat("%s0.dat", pfx ) );
    b->adjust_contour();
    
    
}

YOCTO_UNIT_TEST_IMPL(ctrl)
{
    Real   lam = 0.5;
    //Real   mu  = 0.15;
    Real   gam = 1;
    Bubble bubble(lam,gam,1);
    
    std::cerr << "Circle" << std::endl;
    Shape::Circle(&bubble, Vertex(0,0), 4);
    check(&bubble, "c");
    
    
    
    Shape::Ellipse(&bubble, Vertex(0,0), Vertex(5,2));
    check(&bubble,"e");
    
    Shape::Blob(&bubble, Vertex(0,0), 4, 0.5, 0.0);
    check(&bubble,"b");

}
YOCTO_UNIT_TEST_DONE()
