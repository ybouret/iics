#include "yocto/utest/run.hpp"
#include "../bubbles.hpp"

YOCTO_UNIT_TEST_IMPL(bubble)
{
    const PBC pbc(1);
    Bubbles   bubbles(pbc);
    bubbles.create(10);
    
    for( Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next )
    {
        const Vertex center( 10*(0.5 - Alea()), 10*(0.5-Alea()));
        const Real   radius = 0.1 + 10 * Alea();
        bubble->map_circle(center, radius);
    }
    
}
YOCTO_UNIT_TEST_DONE()
