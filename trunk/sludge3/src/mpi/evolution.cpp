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
    save_markers(MPI);
    for( Bubble *b = bubbles.head; b; b=b->next )
    {
        ios::ocstream fp( vformat("v%u.dat", b->UID), false);
        
        for( Marker *m= b->markers.head; m; m=m->next )
        {
            Tracer       *tr = m->tracer;
            // recompose the gradient
            const Vertex  g  = m->gt * tr->t + m->gn * tr->n;
            
            fp("%g %g\n", tr->pos.x, tr->pos.y);
            fp("%g %g\n\n", tr->pos.x+g.x, tr->pos.y+g.y);
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
