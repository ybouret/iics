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
    
}

void Rescaler:: build_metrics( Bubble &bubble )
{
    const size_t n = bubble.size;
    s.make(n,0);
    ax.make(n,0);
    ay.make(n,0);
    theta.make(n,0);
    Tracer *p = bubble.root;
    period    = 0;
    for( size_t i=1; i <=n; ++i, p=p->next )
    {
        bubble.pbc(p->vertex);
        s[i]  = period;
        ax[i] = p->vertex.x;
        ay[i] = p->vertex.y;
        Vertex pq(p->vertex,p->next->vertex);
        bubble.pbc(pq);
        p->edge = pq;
        p->s2   = pq.norm2();
        p->s    = Sqrt(p->s2);
        period += p->s;
    }
    //! corresponding are
    bubble.compute_area();
    
    //! keep the pressure
    bubble.content = bubble.pressure * bubble.area;
    
}

void Rescaler::abscissa::reset() throw()
{
    s = 0;
    next = prev = 0;
}


void Rescaler:: rebuild( Bubble &bubble )
{
    const size_t n = bubble.size;
    assert( s.size()  == n );
    assert( ax.size() == n );
    assert( ay.size() == n );
    assert( theta.size() == n );
    assert(period>0);
    assert(bubble.area>0);
    assert(a_list.size >= 3);
    
    
    //--------------------------------------------------------------------------
    // build the trigonometric interpolator
    //--------------------------------------------------------------------------
    const Real tfac     = numeric<Real>::two_pi / period;
    for( size_t i=n;i>0;--i)
        theta[i] = tfac * s[i];
    trigonometric<Real> trig( theta, solver);
    
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
    
    //--------------------------------------------------------------------------
    // rebuild its metrics
    //--------------------------------------------------------------------------
    build_metrics(bubble);
    
}
