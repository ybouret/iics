#include "workspace.hpp"

static 
const Junction * __find_next_horz( const Junction::List &JL, const Real x)
{
    assert(JL.size>0);
    bool  in_bubble   = true;
    const Junction *J = JL.head;
    const Junction *K = J->next;
    while(K)
    {
        if(J->value>x)
            return J;
        in_bubble = !in_bubble;
        J=K;
        K=K->next;
    }
    throw exception("missing horizontal junction");
}

void Workspace:: pressurize_contours()
{
    Enter.ldz();
    Leave.ldz();
    //==========================================================================
    // using horizontal axis => x components of Leave/Enter
    //==========================================================================
    for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
    {
        const Junction::List &JL = junctions.Horz(j);
        
        if(JL.size>0)
        {
            
            //------------------------------------------------------------------
            // special case: left side
            //------------------------------------------------------------------
            assert(B[j][lower.x]<0);
            if(B[j][bulk_imin]>=0)
            {
                Enter[j][bulk_imin].x = 0;
            }
            
            //------------------------------------------------------------------
            // bulk case
            //------------------------------------------------------------------
            for( unit_t i=bulk_imin;i<=bulk_imax;++i)
            {
                if(B[j][i]<0)
                {
                    const Real      Xi    = X[i];
                    const Junction *Jnext = B[j][i+1] >= 0 ? __find_next_horz(JL, Xi) : 0;
                    if(Jnext)
                    {
                        
                    }
                }
            }
            
            
            
        }
        
    }
}

#if 0
void Workspace:: pressurize_contours()
{
    //Enter.ldz();
    //Leave.ldz();
    
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
                            const Real phi = J->value - X[J->lower];
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
#endif

