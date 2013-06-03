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
    
    MPI.Printf0(stderr, "\t...validation\n");
    validate_bubbles(MPI);
    if( !is_valid )
    {
        MPI.Printf0(stderr, "Invalid bubble\n");
        done = true;
        return;
    }
    
    MPI.Printf0(stderr, "\t...broadcast bubbles\n");
    broadcast_bubbles(MPI);
    
    MPI.Printf(stderr," #bubbles= %u\n", unsigned(bubbles.size) );
    for(const Bubble *b = bubbles.head;b;b=b->next)
    {
        MPI.Printf(stderr, "\t#%3u = %u\n", unsigned(b->UID), unsigned(b->size) );
    }
    MPI.Printf0(stderr, "\t...segmentation\n");
    segment();
    
    MPI.Printf0(stderr, "\t...pressurize\n");
    P.ldz();
    pressurize_bubbles();
    pressurize_contours();
    compute_gradP(MPI);
    
    
}
