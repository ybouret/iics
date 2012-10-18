#include "../simulation.hpp"

void Simulation:: step()
{
    VisIt::Simulation::step();
    
    //-- advect spots
#if 0
    for( Bubble *b = bubbles.first(); b; b=b->next)
    {
        const Real dv = -0.05 * b->id;
        const Real da = 0.05 * b->id;
        b->rotate(da);

        for( Spot *s = b->spots.head;s;s=s->next)
        {
            Tracer *tracer = s->handle;
            Vertex &v      = tracer->vertex;
            
            
            v.y -= dv;
        }
    }
#else
    for( Bubble *b = bubbles.first(); b; b=b->next)
    {
        for( Spot *s = b->spots.head;s;s=s->next)
        {
            Tracer      *tracer = s->handle;
            Vertex       &v      = tracer->vertex;
            const Vertex &u      = s->U;
            
            v += 0.01 * u;
        }

    }
#endif
    
    
    //-- send info to master
    bubbles.assemble(MPI);

    // TODO: remove this
    for( Bubble *b = bubbles.first(); b; b=b->next)
    {
        //b->set_pressure(1);
    }

    
    //-- legalize new config
    legalize(MPI);
        
}
