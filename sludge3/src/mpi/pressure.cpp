#include "workspace.hpp"


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
        for(unit_t i=outline.lower.x;i<=outline.upper.x;++i)
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
    //Enter.ldz();
    //Leave.ldz();
    return;
    
    //==========================================================================
    // using horizontal axis => x components of Leave/Enter
    //==========================================================================
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
                    if(J->inside||K->inside)
                    {
                        const unit_t ini = J->inside ? J->upper : B.lower.x;
                        const unit_t end = K->inside ? K->lower : B.upper.x;
                        if(end>=ini)
                        {
                            Enter[j][ini].x = J->pressure;
                            Leave[j][end].x = K->pressure;
                        }
                    }
                    
                }
                in_bubble = !in_bubble;
                J=K;
                K=K->next;
            }
        }
    }
    
    
    //==========================================================================
    // using vertical axis => y components of Leave/Enter
    //==========================================================================
    for(unit_t i=outline.lower.x;i<=outline.upper.x;++i)
    {
        const Junction::List &JL = junctions.Vert(i);
        
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
                    
                    if(J->inside||K->inside)
                    {
                        const unit_t ini = J->inside ? J->upper : B.lower.y;
                        const unit_t end = K->inside ? K->lower : B.upper.y;
                        if(end>=ini)
                        {
                            Enter[ini][i].y = J->pressure;
                            Leave[end][i].y = K->pressure;
                        }
                    }
                    
                }
                in_bubble = !in_bubble;
                J=K;
                K=K->next;
            }
        }
    }
    
    
    
    
    
}
