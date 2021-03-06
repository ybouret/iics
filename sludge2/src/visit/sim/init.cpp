#include "../simulation.hpp"

void Simulation:: initialize()
{
    
    bubbles.clear();
    if( is_first )
    {
        //-- create bubbles
        {
            const Vertex center(full_length.x/4,0);
            const Real   radius(full_length.x/9);
            Bubble *bubble = bubbles.append();
            bubble->map_peanut(center,radius,0.9);
            bubble->set_pressure(1.0);
        }
        
        if(true)
        {
            const Vertex center(3*full_length.x/4,full_length.y/4);
            const Real   radius(full_length.x/9);
            Bubble *bubble = bubbles.append();
            bubble->map_peanut(center, radius,0.85);
            bubble->set_pressure(1.0);
            bubble->rotate(M_PI/4);
        }
        
#if 0
        //----------------------------------------------------------------------
        //-- arrange bubbles
        //----------------------------------------------------------------------
        bubbles.update_topology();
#endif
    }
    
    legalize(MPI);
    
}
