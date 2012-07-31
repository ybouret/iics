#include "rescaler.hpp"

Rescaler:: ~Rescaler() throw()
{
}

Rescaler:: Rescaler() :
solver(),
s(),
ax(),
ay(),
period(0),
theta(),
a_pool(),
a_list(a_pool)
{
}

void Rescaler:: process(Bubbles &bubbles)
{
    
    //--------------------------------------------------------------------------
    // phase 1: rescale every bubbles to match their critical parameters
    //--------------------------------------------------------------------------
    for( Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next )
    {
        upgrade(*bubble);
    }
    
}

void Rescaler:: build_metrics( Bubble &bubble )
{
    const size_t n = bubble.size;
    std::cerr << "building metrics for " << n << " points" << std::endl;
    s.make(n,0);
    ax.make(n,0);
    ay.make(n,0);
    theta.make(n,0);
    Tracer *p = bubble.root;
    period    = 0;
    Vertex  v0 = p->vertex;
    Real    area = 0;
    for( size_t i=1; i <=n; ++i, p=p->next )
    {
        bubble.pbc(p->vertex);
        s[i]  = period;
        
        Vertex pq(p->vertex,p->next->vertex);
        bubble.pbc(pq);
        p->edge = pq;
        p->s2   = pq.norm2();
        p->s    = Sqrt(p->s2);
        period += p->s;
        ax[i] = v0.x;
        ay[i] = v0.y;
        const Vertex v1 = v0 + pq;
        area += v0.x * v1.y - v0.y * v1.x;
        v0 = v1;
    }
    //! corresponding area
    bubble.area = Fabs(area)/2;
    std::cerr << "\tarea=" << bubble.area << std::endl;
    
    //! keep the pressure
    bubble.content = bubble.pressure * bubble.area;
    
}

void Rescaler::abscissa::reset() throw()
{
    s    = 0;
    next = prev = 0;
}

#include "yocto/ios/ocstream.hpp"

void Rescaler:: rebuild( Bubble &bubble )
{
    static int fid = 0;
    
    const size_t n = bubble.size;
    assert( s.size()     == n );
    assert( ax.size()    == n );
    assert( ay.size()    == n );
    assert( theta.size() == n );
    assert(period>0);
    assert(bubble.area>0);
    assert(a_list.size >= 3);
    
    
    {
        ios::ocstream fp( vformat("org%d.dat",fid), false);
        for( size_t i=1; i <= n; ++i )
        {
            fp("%g %g %g\n", s[i], ax[i], ay[i]);
        }
        fp("%g %g %g\n", period, ax[1], ay[1]);
    }
    
    //--------------------------------------------------------------------------
    // build the trigonometric interpolator
    //--------------------------------------------------------------------------
    const Real tfac     = numeric<Real>::two_pi / period;
    for( size_t i=n;i>0;--i)
        theta[i] = tfac * s[i];
    trigonometric<Real> trig( theta, solver);
    trig.compute(ax, solver);
    trig.compute(ay, solver);
    
    //--------------------------------------------------------------------------
    // rebuild the bubble
    //--------------------------------------------------------------------------
    bubble.empty();
    for( const abscissa *a = a_list.head; a; a=a->next )
    {
        const Real s_i = a->s;
        assert(s_i>=0);
        assert(s_i<period);
        const Real arg = s_i * tfac;
        bubble.append()->vertex = trig( arg, ax, ay );
    }
    assert(bubble.size==a_list.size);
    
    {
        ios::ocstream fp( vformat("ref%d.dat",fid), false);
        const abscissa *a = a_list.head;
        const Tracer   *p = bubble.root;
        while(a)
        {
            fp("%g %g %g\n", a->s, p->vertex.x, p->vertex.y);
            a=a->next;
            p=p->next;
        }
        fp("%g %g %g\n", period, p->vertex.x, p->vertex.y);

        ++fid;
    }
    
    //--------------------------------------------------------------------------
    // rebuild its metrics
    //--------------------------------------------------------------------------
    build_metrics(bubble);
    
}

void Rescaler:: upgrade( Bubble &bubble )
{
    build_metrics(bubble);
    refine(bubble);
}
