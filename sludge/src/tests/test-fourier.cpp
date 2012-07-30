#include "yocto/utest/run.hpp"
#include "../bubble.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/kernel/linsys.hpp"
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



static inline
void scale_fourier( Bubble &bubble )
{
    vector<Real> s;
    vector<Real> a;
    vector<Real> b;
    
    
    const size_t    n = bubble.size;
    const bool      is_even = 0 == (n&1);
    const size_t    nn = (n >> 1) - ( is_even ? 1 : 0 );
    std::cerr << "n=" << n << ", nn=" << nn << std::endl;
    Tracer *p      = bubble.root;
    Real    period = 0;
    for( size_t i=n;i>0;--i,p=p->next)
    {
        bubble.pbc(p->vertex);
        {
            s.push_back(period);
            
            a.push_back(p->vertex.x);
            b.push_back(p->vertex.y);
            
        }
        Vertex edge(p->vertex,p->next->vertex);
        bubble.pbc(edge);
        period += edge.norm();
    }
    
    const double tfac = numeric<Real>::two_pi/period;
    
    
    matrix<Real> M(n,n);
    for( size_t i=1; i <=n; ++i )
    {
        const Real t_i = s[i] * tfac;
        size_t j       = 1;
        M[i][j++] = 1;
        for( size_t k=1; k <=nn; ++k )
        {
            const Real    arg = k * t_i;
            M[i][j++] = Cos(arg);
            M[i][j++] = Sin(arg);
        }
        if(is_even)
            M[i][j] = Cos( (nn+1)*t_i );
    }
    
    {
        ios::ocstream fp("source.dat", false);
        for( size_t i=1; i <= n; ++i )
        {
            fp("%g %g %g\n", s[i]*tfac, a[i], b[i] );
        }
    }
    
    linsys<Real> solver(n);
    if( !solver.LU(M) )
    {
        throw exception("Matrix is not invertible");
    }
    solver(M,a);
    solver(M,b);
    //std::cerr << "a=" << a << std::endl;
    //std::cerr << "b=" << b << std::endl;
    
    const size_t m = max_of<size_t>(3,ceil( period/ bubble.lambda ));
    std::cerr << "Placing " << m << " points from " << n << std::endl;
    const double ds = period/m;
    bubble.empty();
    ios::ocstream fp("target.dat",false);
    const Real t0 = Alea() * numeric<Real>::two_pi;
    for( size_t i=0; i < m; ++i )
    {
        const Real t_i = t0 + (i*ds) * tfac;
        size_t     j=1;
        
        Real x = a[j];
        Real y = b[j];
        ++j;
        for( size_t k=1; k <= nn; ++k )
        {
            const Real    arg = k * t_i;
            const Real    ca  = Cos(arg);
            const Real    sa  = Sin(arg);
            x += a[j]   * ca;
            y += b[j]   * ca;
            ++j;
            x += a[j]   * sa;
            y += b[j]   * sa;
            ++j;
        }
         if(is_even)
         {
             const Real arg = (nn+1) * t_i;
             x += a[j] * Cos( arg );
             y += b[j] * Sin( arg );
         }
        fp("%g %g %g\n", t_i, x, y );
        Vertex &v = bubble.append()->vertex;
        v.x = x;
        v.y = y;
    
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
    
    bubble.map_circle( center, 2.0);
    bubble.save_dat("circle0.dat");
    scale_fourier(bubble);
    bubble.save_dat("circle1.dat");
    expand(bubble,center);
    bubble.save_dat("circle2.dat");
    scale_fourier(bubble);
    bubble.save_dat("circle3.dat");
    
    bubble.map_peanut( center, 2, 0.96);
    bubble.save_dat("peanut0.dat");
    scale_fourier(bubble);
    bubble.save_dat("peanut1.dat");

    expand(bubble,center);
    bubble.save_dat("peanut2.dat");
    scale_fourier(bubble);
    bubble.save_dat("peanut3.dat");
    
}
YOCTO_UNIT_TEST_DONE()
