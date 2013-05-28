#include "yocto/utest/run.hpp"
#include "yocto/code/rand.hpp"
#include "shape.hpp"


YOCTO_UNIT_TEST_IMPL(bubble)
{
    Real   lam = 1;
    Bubble bubble(lam);
    const size_t N = 3+alea_leq(50);
    for( size_t i=1; i <= N; ++i )
    {
        Tracer *tr = new Tracer();
        bubble.push_back(  tr );
        const Real theta = ((i-1) * numeric<Real>::two_pi)/(N);
        const Real rho   = Sqrt(5.0+4.0*Cos(3.0*theta));
        tr->pos.x = rho * Cos(theta);
        tr->pos.y = rho * Sin(theta);
    }
    Hasher h;
    h.set();
    bubble.hash_bubble(h);
    const Hasher::KeyType k = h.getKey();
    std::cerr << "key=" << k << std::endl;
    
    bubble.save_dat( "bubble.dat" );
    bubble.init_contour();
    bubble.auto_contour();
    bubble.save_dat( "autob.dat" );
    
    Shape::Ellipse(&bubble, Vertex(0,0), Vertex(6,2));
    Shape::Rotate(&bubble, 0.2);
    bubble.save_dat("ell.dat");
    bubble.auto_contour();
    bubble.save_dat("ell2.dat");
    
}
YOCTO_UNIT_TEST_DONE()
