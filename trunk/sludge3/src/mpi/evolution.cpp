#include "workspace.hpp"


void Workspace:: evolution(const mpi &MPI, Real dt)
{
    
    //==========================================================================
    //
    // Evolution of all the data fields
    //
    //==========================================================================
    
    
    
    
    //==========================================================================
    //
    // Evolution of the markers
    //
    //==========================================================================
    MPI.Printf0(stderr, "\t\tevolving markers...\n");
    for( Bubble *b = bubbles.head; b; b=b->next )
    {
        for( Marker *m= b->markers.head; m; m=m->next )
        {
            Tracer       *tr = m->tracer;
            // recompose the gradient
            const Vertex  g  = m->gt * tr->t + m->gn * tr->n;
            
            // compute the velocity
            const Vertex  v  = gradP_to_V(g);
            
            // deduce the displacement
            const Vertex  dr = dt * v;

            // evolve
            tr->pos += dr;
        }
    }
    
    //==========================================================================
    //
    // Send back markers to the master
    //
    //==========================================================================
    MPI.Printf0(stderr, "\t\tcollecting markers....\n");
    ParallelBubbles::Collect(MPI, bubbles);
    
}
