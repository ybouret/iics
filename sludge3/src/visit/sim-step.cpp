#include "simulation.hpp"
#include "../shape.hpp"

void Simulation:: step()
{
    
    MPI.Printf0(stderr, "[SIMULATION STEP]\n");
#if 0
    if(MPI.IsFirst)
    {
        for( Bubble *bubble = bubbles.head;bubble;bubble=bubble->next)
        {
            const Real omega = 0.02 * (bubble->UID+1);
            Shape::Rotate(bubble, omega);
        }
    }
#endif
    
    evolution(MPI, 0.01);
    initialize();
    
    
}
