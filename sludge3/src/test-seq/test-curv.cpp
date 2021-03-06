#include "yocto/utest/run.hpp"
#include "yocto/code/rand.hpp"
#include "shape.hpp"


YOCTO_UNIT_TEST_IMPL(curv)
{
    Real   lam = 1;
    Real   gam = 1;
    Bubble bubble(lam,gam,2);
    
    
    Shape::Ellipse(&bubble, Vertex(0,0), Vertex(6,2));
    Shape::Rotate(&bubble, 0.3);
    
    bubble.adjust_contour();
    bubble.compute_curvatures();
    
    bubble.save_dat("ell.dat");
    bubble.save_t("ell_t.dat");
    bubble.save_n("ell_n.dat");
    
    
    Shape::Circle(&bubble, Vertex(0,0), 5);
    bubble.adjust_contour();
    bubble.compute_curvatures();
    
    bubble.save_all("circ");    
    
    std::cerr << "sizeof(Tracer)=" << sizeof(Tracer) << std::endl;

}
YOCTO_UNIT_TEST_DONE()
