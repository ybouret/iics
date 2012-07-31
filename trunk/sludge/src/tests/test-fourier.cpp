#include "yocto/utest/run.hpp"
#include "../bubble.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"

static inline
void expand( Bubble &bubble , const Vertex &center)
{
    Tracer *p = bubble.root;
    for( size_t i=bubble.size;i>0;--i,p=p->next)
    {
        const Real    theta = p->vertex.angle();
        const Vertex  v( 1.5 * Cos(theta), 2.0 * Cos( theta ) );
        const Vertex  r = p->vertex - center;
        p->vertex = center + v.norm() * r;
    }
}



#include "yocto/math/dat/trigonometric.hpp"

void scale_trigo( Bubble &bubble )
{
    vector<Real>    s;
    vector<Real>    ax;
    vector<Real>    ay;
    
    const size_t n      = bubble.size;
    Tracer      *p      = bubble.root;
    Real         period = 0;
    std::cerr << "trigo for " << n << " points" << std::endl;
    
    for( size_t i=n;i>0;--i,p=p->next)
    {
        bubble.pbc(p->vertex);
        {
            s.push_back(period);
            
            ax.push_back(p->vertex.x);
            ay.push_back(p->vertex.y);
            
        }
        Vertex edge(p->vertex,p->next->vertex);
        bubble.pbc(edge);
        period += edge.norm();
    }

    const Real   tfac = numeric<Real>::two_pi / period;
    vector<Real> theta(s);
    for( size_t i=n;i>0;--i) theta[i] *= tfac;
    linsys<Real>        solver;
    trigonometric<Real> trig( theta, solver );
    trig.compute(ax, solver);
    trig.compute(ay, solver);
    
    
    const size_t m = 5*max_of<size_t>(3,ceil( period/ bubble.lambda ));
    std::cerr << "Placing " << m << " points from " << n << std::endl;
    const double ds = period/m;
    bubble.empty();

    for( size_t i=0; i < m; ++i )
    {
        const Real t_i = (i*ds) * tfac;
        bubble.append()->vertex = trig( t_i, ax, ay);
    }
    
}

YOCTO_UNIT_TEST_IMPL(fourier)
{
    AleaInit();
    double lambda = 1;
    Vertex box(100,100);
    PBC    pbc(box.y);
    Tracer::Cache tcache;
    Spot::Cache   scache;
    Marker::Cache mcache;
    
    Vertex center( box.x/2, 0.0 );
    
    Bubble bubble(lambda,pbc,tcache,scache,mcache);
    
    bubble.map_circle( center, 2.0 + Alea() * 0.2 );
    bubble.save_dat("circle0.dat");
    scale_trigo(bubble);
    bubble.save_dat("circle1.dat");
    expand(bubble,center);
    bubble.save_dat("circle2.dat");
    scale_trigo(bubble);
    bubble.save_dat("circle3.dat");
    
    bubble.map_peanut( center, 2+Alea()*0.2, 0.96);
    bubble.save_dat("peanut0.dat");
    scale_trigo(bubble);
    bubble.save_dat("peanut1.dat");
    expand(bubble,center);
    bubble.save_dat("peanut2.dat");
    scale_trigo(bubble);
    bubble.save_dat("peanut3.dat");
}
YOCTO_UNIT_TEST_DONE()
