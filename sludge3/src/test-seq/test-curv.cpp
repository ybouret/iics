#include "yocto/utest/run.hpp"
#include "yocto/code/rand.hpp"
#include "shape.hpp"


YOCTO_UNIT_TEST_IMPL(curv)
{
    Real   lam = 1;
    Bubble bubble(lam);
    
    
    Shape::Ellipse(&bubble, Vertex(0,0), Vertex(6,2));
    Shape::Rotate(&bubble, 0.3);
    
    bubble.auto_contour();
    bubble.compute_curvatures();
    
    bubble.save_dat("ell.dat");
    bubble.save_t("ell_t.dat");
    bubble.save_n("ell_n.dat");
    
    
    

}
YOCTO_UNIT_TEST_DONE()
