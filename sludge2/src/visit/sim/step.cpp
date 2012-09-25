#include "../simulation.hpp"

void Simulation:: step()
{
    VisIt::Simulation::step();
    
    //-- advect spots
    for( Bubble *b = bubbles.first(); b; b=b->next)
    {
        for( Spot *s = b->spots.head;s;s=s->next)
        {
            Tracer *tracer = s->handle;
            Vertex &v      = tracer->vertex;
            v.y += 0.1;
        }
    }
    
    
    //-- send info to master
    bubbles.assemble(MPI);
    
    //-- legalize new config
    legalize(MPI);
    
    
    
}