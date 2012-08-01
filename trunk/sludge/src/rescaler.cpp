#include "rescaler.hpp"

Rescaler:: ~Rescaler() throw()
{
}

Rescaler:: Rescaler() :
s(),
ax(),
ay(),
period(0),
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
    const size_t n1 = n+1;
    std::cerr << "building metrics for " << n << " points" << std::endl;
    s.make(n1,0);
    ax.make(n1,0);
    ay.make(n1,0);
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
    s[n1]  = period;
    ax[n1] = ax[1];
    ay[n1] = ay[1];
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


void Rescaler:: upgrade( Bubble &bubble )
{
    build_metrics(bubble);
    refine(bubble);
}
