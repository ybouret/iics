#include "workspace.hpp"

Vertex Workspace:: gradP_to_V( const Vertex &g ) const
{
    return -g;
}


void Workspace:: compute_velocities()
{
    
    //==========================================================================
    //
    // compute velocities in the bulk
    //
    //==========================================================================
    for(unit_t j=outline.lower.y;j<=outline.upper.y;++j)
    {
        for(unit_t i=outline.lower.x;i<=outline.upper.x;++i)
        {
            if(B[j][i] < 0 )
            {
                V[j][i] = gradP_to_V(gradP[j][i]);
            }
            else
                V[j][i].ldz();
        }
        
    }
    
    //==========================================================================
    //
    // compute for all the markers of all the bubbles
    //
    //==========================================================================
    for( Bubble *b = bubbles.head;b;b=b->next)
    {
        const Real P_in  = b->pressure;
        const Real gamma = b->gamma;
        for( Marker *m = b->markers.head;m;m=m->next)
        {
            //------------------------------------------------------------------
            // Tangential, easy
            //------------------------------------------------------------------

            const Tracer *tr    = m->tracer;
            const Real    Pcurr = P_in - gamma * tr->C;
            const Real    Pnext = P_in - gamma * tr->next->C;
            const Real    Pprev = P_in - gamma * tr->prev->C;
            const Real    Pp    = Pnext-Pcurr;
            const Real    Pm    = Pprev-Pcurr;
            const Real    tp    = tr->dist;
            const Real    tm    = tr->prev->dist;
            const Real    Vp    = Pp/tp;
            const Real    Vm    = Pm/tm;
            m->gt = (tm*Vp - tp*Vm)/(tm+tp);
            
            //------------------------------------------------------------------
            // normal pressure
            //------------------------------------------------------------------
            m->gn = 0;
        }
    }
    
    
}



