#include "yocto/utest/run.hpp"
#include "yocto/code/rand.hpp"
#include "shape.hpp"


YOCTO_UNIT_TEST_IMPL(curv)
{
    Real   lam = 1;
    Bubble bubble(lam,2);
    
    
    Shape::Ellipse(&bubble, Vertex(0,0), Vertex(6,2));
    Shape::Rotate(&bubble, 0.3);
    
    bubble.auto_contour();
    bubble.compute_curvatures();
    
    bubble.save_dat("ell.dat");
    bubble.save_t("ell_t.dat");
    bubble.save_n("ell_n.dat");
    
    
    Shape::Circle(&bubble, Vertex(0,0), 5);
    bubble.auto_contour();
    bubble.compute_curvatures();
    
    bubble.save_dat("circ.dat");
    bubble.save_t("circ_t.dat");
    bubble.save_n("circ_n.dat");
    
    
    std::cerr << "sizeof(Tracer)=" << sizeof(Tracer) << std::endl;

}
YOCTO_UNIT_TEST_DONE()
