#include "workspace.hpp"
#include "yocto/code/utils.hpp"

bool Workspace:: evolution(const mpi &MPI, Real dt)
{
    
    //==========================================================================
    //
    // assume the velocities/pressure etc.. are computed
    //
    //==========================================================================
    save_markers(MPI);
    while( dt > 0 )
    {
        Real dt_max = dt;
        
        //----------------------------------------------------------------------
        // compute the velocities for each marker, and
        // deduced max allowed time
        //----------------------------------------------------------------------
        for( Bubble *b = bubbles.head; b; b=b->next )
        {
            //ios::ocstream fp( vformat("v%u.dat", unsigned(b->UID) ), false);
            
            for( Marker *m= b->markers.head; m; m=m->next )
            {
                Tracer       *tr = m->tracer;
                // recompose the gradient
                const Vertex  g  = m->gt * tr->t + m->gn * tr->n;
                
                //fp("%g %g\n", tr->pos.x, tr->pos.y);
                //fp("%g %g\n\n", tr->pos.x+g.x, tr->pos.y+g.y);
                
                // compute the velocity
                m->v  = gradP_to_V(g);
                
                // compute the max allowed time
                const Real d_max = min_of<Real>( tr->dist, tr->prev->dist)/3;
                const Real speed = m->v.norm();
                if( speed * dt_max > d_max )
                {
                    dt_max = d_max / speed;
                }
            }
            
            
        }
        
        //----------------------------------------------------------------------
        // move according to the min of dt_max
        //----------------------------------------------------------------------
        MPI.Printf(stderr,"dt_max=%g\n", dt_max);
        dt_max = MPI.Min(dt_max,MPI_COMM_WORLD);
        MPI.Printf0(stderr, "\tdt_max=%g\n", dt_max);
        
        for( Bubble *b = bubbles.head; b; b=b->next )
        {
            for( Marker *m= b->markers.head; m; m=m->next )
            {
                Tracer  *tr = m->tracer;
                tr->pos += dt_max * m->v;
            }
        }
        dt -= dt_max;
        
        MPI.Printf0(stderr, "\t\tcollecting markers....\n");
        ParallelBubbles::Collect(MPI, bubbles);
        
        if(!initialize(MPI))
            return false;
    }
    
    
    return true;
}
