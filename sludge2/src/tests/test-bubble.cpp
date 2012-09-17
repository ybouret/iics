#include "yocto/utest/run.hpp"
#include "../bubbles.hpp"

YOCTO_UNIT_TEST_IMPL(bubble)
{
    const PBC pbc(50);
    Bubbles   bubbles(pbc);
    bubbles.create(10);
    
    for( Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next )
    {
        const Vertex center( 10*(0.5 - Alea()), 10*(0.5-Alea()));
        const Real   radius = 0.1 + 10 * Alea();
        bubble->map_peanut(center, radius,0.9);
        if( bubble == bubbles.first())
        {
            std::cerr << "#tracers=" << bubble->size << std::endl;
            bubble->save_dat( "b0.dat" );
            bubble->compute_contour();
            bubble->save_dat( "b1.dat" );
            bubble->save_vtk(   "b1.vtk");
            bubble->save_vtk_t( "b1_t.vtk" );
            bubble->save_vtk_n( "b1_n.vtk" );
        }
    }
    
}
YOCTO_UNIT_TEST_DONE()
