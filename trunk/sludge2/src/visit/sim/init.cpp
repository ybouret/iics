#include "../simulation.hpp"

void Simulation:: initialize()
{
    
    bubbles.clear();
    if( master )
    {
        //-- create bubbles
        const Vertex center(full_length.x/3,0);
        const Real   radius(full_length.x/6);
        Bubble *bubble = bubbles.append();
        bubble->map_peanut(center,radius,0.9);
        bubble->set_pressure(1.0);
        
        //----------------------------------------------------------------------
        //-- arrange bubbles
        //----------------------------------------------------------------------
        bubbles.update_topology();
    }
    
    legalize(MPI);
    
}
