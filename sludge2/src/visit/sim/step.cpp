#include "../simulation.hpp"

void Simulation:: step()
{
    VisIt::Simulation::step();
    
    //--------------------------------------------------------------------------
    //-- advect spots
    //--------------------------------------------------------------------------
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
    MPI.Printf0(stderr,"\tadvecting spots/tracers\n");
    for( Bubble *b = bubbles.first(); b; b=b->next)
    {
        for( Spot *s = b->spots.head;s;s=s->next)
        {
            Tracer      *tracer = s->handle;
            Vertex       &v     = tracer->vertex;
            const Vertex &u     = s->U;
            //fprintf(stderr,"(%g,%g) => ", v.x,v.y);
            v += 0.02 * u;
            //fprintf(stderr,"(%g,%g)\n", v.x,v.y);
        }
    }
#endif
    
    //--------------------------------------------------------------------------
    //-- send info to master
    //--------------------------------------------------------------------------
    MPI.Printf0(stderr,"\tassembling bubbles\n");
    bubbles.assemble(MPI);
    
    
    // TODO: remove this
    for( Bubble *b = bubbles.first(); b; b=b->next)
    {
        //b->set_pressure(1);
    }
    
    
    //-- legalize new config
    legalize(MPI);
    
}
