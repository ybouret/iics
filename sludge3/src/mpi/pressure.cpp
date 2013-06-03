#include "workspace.hpp"

void Workspace:: reset_pressure()
{
    P.ldz();
    pressurize_bubbles();
}

void Workspace:: pressurize_bubbles()
{
    //--------------------------------------------------------------------------
    // Collect the pressures
    //--------------------------------------------------------------------------
    bpres.free();
    bpres.reserve(bubbles.size);
    for( const Bubble *b = bubbles.head;b;b=b->next)
    {
        bpres.push_back( b->pressure );
    }
    
    //--------------------------------------------------------------------------
    // Fill the bubbles location with bubble pressures
    //--------------------------------------------------------------------------
    for( unit_t j=outline.lower.y; j <= outline.upper.y; ++j)
    {
        for(unit_t i=lower.x;i<=upper.x;++i)
        {
            const Real which = B[j][i];
            if( which >= 0 )
            {
                assert(which<bubbles.size);
                assert(which<bpres.size());
                P[j][i] = bpres[ which+1 ];
            }
            Enter[j][i].x = Enter[j][i].y =
            Leave[j][i].x = Leave[j][i].y = P[j][i];
        }
    }
    
    
}


void Workspace:: pressurize_contours()
{
    //! using horizontal axis => x components of Leave/Enter
    for(unit_t j=outline.lower.y;j<=outline.upper.y;++j)
    {
        const Junction::List &JL = junctions.Horz(j);
        
        if(JL.size>0)
        {
            bool in_bubble = true;
            const Junction *J = JL.head;
            const Junction *K = J->next;
            while(K)
            {
                if(in_bubble)
                {
                    if(J->owner!=K->owner)
                        throw exception("Bubble in Bubble!");
                    
                }
                in_bubble = !in_bubble;
                J=K;
                K=K->next;
            }
        }
    }
    
    
    
}
