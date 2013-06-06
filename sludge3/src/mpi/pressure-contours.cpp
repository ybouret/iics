#include "workspace.hpp"


void Workspace:: EnterX(const Junction *J, unit_t j)
{
    assert(J);
    assert(J->active);
    assert(Bubble::IsAfter == J->b_pos);
    
    const unit_t i = J->lower;
    assert(B[j][J->upper]>=0); // in bubble since J->active
    if(i>=bulk_imin)
    {
        const Real Pm  = P[j][i-1];
        const Real Xi  = X[i];
        const Real phi = J->value - Xi; assert(phi>=0);
        const Real Pb  = J->pressure;
        const Real sigma_E = (Pb - Pm) / (delta.x + phi);
        Enter[j][J->upper].x = Pm + two_delta.x * sigma_E;
    }
}


void Workspace:: LeaveX(const Junction *K, unit_t j)
{
    assert(K);
    assert(K->active);
    assert(Bubble::IsBefore == K->b_pos);
    
    assert(B[j][K->lower]>=0); // since K->active
    const unit_t i = K->upper;
    if(i<=bulk_imax)
    {
        const Real Pp    = P[j][i+1];
        const Real Xi    = X[i];
        const Real psi   = Xi - K->value; assert(psi>=0);
        const Real Pb    = K->pressure;
        const Real sigma_L = (Pp - Pb) / (delta.x + psi);
        Leave[j][K->lower].x = Pp - two_delta.x *sigma_L;
    }
}


void Workspace:: AloneX(const Junction *J, const Junction *K, unit_t j)
{
    assert(J);
    assert(Bubble::IsBefore == J->b_pos);
    assert(J->active);
    assert(B[j][J->lower]>=0); // since J->active
    
    assert(K);
    assert(K==J->next);
    assert(Bubble::IsAfter == K->b_pos);
    assert(K->active);
    assert(B[j][K->upper]>=0); // since K->active
    assert(J->upper == K->lower);
    
    const Real fac= two_delta.x / (K->value - J->value);
    Leave[j][J->lower].x = J->pressure * fac;
    Enter[j][K->upper].x = K->pressure * fac;
    
    
}

void Workspace:: pressurize_horz()
{
    
    for( unit_t j=outline.lower.y; j<=outline.upper.y; ++j)
    {
        const Junction::List &JL = junctions.Horz(j);
        if(JL.size<=1) continue;
        
        //----------------------------------------------------------------------
        // before crossing left-most bubble.
        //----------------------------------------------------------------------
        const Junction *J = JL.head;
        assert(Bubble::IsAfter == J->b_pos);
        if(J->active)
        {
            EnterX(J,j);
        }
        
        //----------------------------------------------------------------------
        // between bubbles
        //----------------------------------------------------------------------
        bool  in_bubble   = true;
        const Junction *K = J->next;
        while(K)
        {
            
            if(!in_bubble)
            {
                assert(Bubble::IsBefore == J->b_pos);
                assert(Bubble::IsAfter  == K->b_pos);
                if(J->active)
                {
                    if(K->active)
                    {
                        const unit_t ini = J->upper;
                        const unit_t end = K->lower;
                        if(end>=ini)
                        {
                            if(end>ini)
                            {
                                LeaveX(J, j);
                                EnterX(K, j);
                            }
                            else
                            {
                                assert(ini==end);
                                AloneX(J,K,j);
                            }
                        }
                        // else pass-through...
                    }
                    else
                    {
                        LeaveX(J, j);
                    }
                }
                else
                {
                    if(K->active)
                    {
                        EnterX(K, j);
                    }
                    else
                    {
                        // do nothing
                    }
                }
            }
            
            in_bubble = !in_bubble;
            J = K;
            K = K->next;
        }
        
        //----------------------------------------------------------------------
        // after crossing last bubble
        //----------------------------------------------------------------------
        assert(!in_bubble);
        assert(J!=0);
        if(J->active)
        {
            LeaveX(J,j);
        }
        
        
    }
    
    
}


void Workspace:: EnterY(const Junction *J, unit_t i)
{
    assert(J);
    assert(J->active);
    assert(Bubble::IsAfter == J->b_pos);
    
    const unit_t j = J->lower;
    assert(B[J->upper][i]>=0); // in bubble since J->active
    if(j>=bulk_jmin)
    {
        const Real Pm  = P[j-1][i];
        const Real Yj  = Y[j];
        const Real phi = J->value - Yj; assert(phi>=0);
        const Real Pb  = J->pressure;
        Enter[J->upper][i].y = Pm + two_delta.y * (Pb - Pm) / (delta.y + phi);
    }
}


void Workspace:: LeaveY(const Junction *K, unit_t i)
{
    assert(K);
    assert(K->active);
    assert(Bubble::IsBefore == K->b_pos);
    assert(B[K->lower][i]>=0); //since K->active

    const unit_t j = K->upper;
    if(j<=bulk_jmax)
    {
        const Real Pp  = P[j+1][i];
        const Real Yj  = Y[j];
        const Real psi = Yj - K->value; assert(psi>=0);
        const Real Pb  = K->pressure;
        Leave[K->lower][i].y = Pp + two_delta.y * (Pb - Pp) / (delta.y + psi );
    }
}


void Workspace:: AloneY(const Junction *J, const Junction *K, unit_t i)
{
    
    assert(J);
    assert(Bubble::IsBefore == J->b_pos);
    assert(J->active);
    
    assert(K);
    assert(K==J->next);
    assert(Bubble::IsAfter == K->b_pos);
    assert(K->active);
    assert(J->upper == K->lower);
    
    assert(B[K->upper][i]>=0); // since K is active
    assert(B[J->lower][i]>=0); // since J is active
    
    const Real fac= two_delta.y / (K->value - J->value);
    Leave[J->lower][i].y = J->pressure * fac;
    Enter[K->upper][i].y = K->pressure * fac;
    
}


void Workspace:: pressurize_vert()
{
    for( unit_t i=outline.lower.x; i<=outline.upper.x; ++i)
    {
        const Junction::List &JL = junctions.Vert(i);
        if(JL.size<=1) continue;
        
        //----------------------------------------------------------------------
        // before crossing left-most bubble.
        //----------------------------------------------------------------------
        const Junction *J = JL.head;
        assert(Bubble::IsAfter == J->b_pos);
        if(J->active)
        {
            EnterY(J,i);
        }
        
        //----------------------------------------------------------------------
        // between bubbles
        //----------------------------------------------------------------------
        bool  in_bubble   = true;
        const Junction *K = J->next;
        while(K)
        {
            
            if(!in_bubble)
            {
                assert(Bubble::IsBefore == J->b_pos);
                assert(Bubble::IsAfter  == K->b_pos);
                if(J->active)
                {
                    if(K->active)
                    {
                        const unit_t ini = J->upper;
                        const unit_t end = K->lower;
                        if(end>=ini)
                        {
                            if(end>ini)
                            {
                                LeaveY(J, i);
                                EnterY(K, i);
                            }
                            else
                            {
                                assert(ini==end);
                                AloneY(J,K,i);
                            }
                        }
                        // else pass-through...
                    }
                    else
                    {
                        LeaveY(J, i);
                    }
                }
                else
                {
                    if(K->active)
                    {
                        EnterY(K, i);
                    }
                    else
                    {
                        // do nothing
                    }
                }
            }
            
            in_bubble = !in_bubble;
            J = K;
            K = K->next;
        }
        
        //----------------------------------------------------------------------
        // after crossing last bubble
        //----------------------------------------------------------------------
        assert(!in_bubble);
        assert(J!=0);
        if(J->active)
        {
            LeaveY(J,i);
        }
        
        
    }

}


void Workspace:: pressurize_contours()
{
    //Enter.ldz();
    //Leave.ldz();
    
    pressurize_horz();
    pressurize_vert();
    
}



