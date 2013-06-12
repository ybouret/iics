#include "simulation.hpp"
#include "../shape.hpp"

void Simulation:: step()
{
    
    MPI.Printf0(stderr, "[SIMULATION STEP]\n");
    MPI.Printf0(stderr, "\t...rotation\n");
    if(MPI.IsFirst)
    {
        for( Bubble *bubble = bubbles.head;bubble;bubble=bubble->next)
        {
            const Real omega = 0.02 * (bubble->UID+1);
            Shape::Rotate(bubble, omega);
        }
    }
    
    validate_bubbles(MPI);
    if(!is_valid)
        throw exception("Invalid Bubbles");
    broadcast_bubbles(MPI);
    segment();
    compute_pressure(MPI, 1e-5);
    
    
    
    
}
