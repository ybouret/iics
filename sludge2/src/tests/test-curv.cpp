#include "yocto/utest/run.hpp"
#include "../bubbles.hpp"


namespace
{
    static inline
    void check( const Bubble &bubble )
    {
        std::cerr << "Bubble #tracers=" << bubble.size << std::endl;
        const Tracer *p = bubble.root;
        for( size_t i=0;i<bubble.size;++i,p=p->next)
        {
            std::cerr << "p=" << p->vertex << "; dsc=" << p->dsc << std::endl;
            std::cerr << "\tt=" << p->t << std::endl;
            std::cerr << "\tC=" << p->curvature << " / "  << 2*Vertex::det(p->t,p->edge) /p->s2 << std::endl;
        }
    }
}

YOCTO_UNIT_TEST_IMPL(curv)
{
    const PBC pbc(50);
    Bubbles   bubbles(pbc);
    bubbles.create(1);
    
    Bubble *bubble = bubbles.first();
    
    const Vertex center(0,0);
    bubbles.lambda = 0.3;
    bubble->map_peanut(center, 1.7, 0.8);
    bubbles.update_topology();
    check( *bubble );
    
    bubble->map_circle(center, 1.7);
    bubbles.update_topology();
    check( *bubble );
    
   

    
}
YOCTO_UNIT_TEST_DONE()
