#include "../simulation.hpp"

void Simulation:: step()
{
    VisIt::Simulation::step();
    
    //-- advect spots
    for( Bubble *b = bubbles.first(); b; b=b->next)
    {
        b->rotate(0.12);
        for( Spot *s = b->spots.head;s;s=s->next)
        {
            Tracer *tracer = s->handle;
            Vertex &v      = tracer->vertex;
            v.y += 0.05;
        }
    }
    
    
    //-- send info to master
    bubbles.assemble(MPI);
    
    //-- legalize new config
    legalize(MPI);
    
    
    
}