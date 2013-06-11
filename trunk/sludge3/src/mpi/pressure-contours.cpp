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
        const unit_t im    = i-1;
        const unit_t ip    = i+1;
        const Real   P_im  = P[j][im];
        const Real   Xi    = X[i];
        const Real   phi   = J->value - Xi; assert(phi>=0);
        const Real   P_phi = J->pressure;
        const Real   L     = delta.x + phi;
        const Real   sigma = (P_phi - P_im) / L;
        E1[j][ip].x        = P_im + two_delta.x * sigma;
        
        
        const Real P_i     = P[j][i];
        const Real h       = L/2;
        const Real dx      = delta.x;
        const Real Peff    = (h*P_i + (dx-h) * P_im)/dx;
        const Real d2P     = P_im - (Peff+Peff) + P_phi;
        const Real fac     = dx/h;
        E2[j][ip].x        = (P_i+P_i) - P_im + fac*fac * d2P;
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
        const unit_t ip    = i+1;
        const unit_t im    = i-1;
        const Real   P_ip  = P[j][ip];
        const Real   Xi    = X[i];
        const Real   psi   = Xi - K->value; assert(psi>=0);
        const Real   P_psi = K->pressure;
        const Real   L     = (delta.x+psi);
        const Real   sigma = (P_ip - P_psi) / (delta.x + psi);
        L1[j][im].x        = P_ip - two_delta.x *sigma;
        
        const Real  P_i    = P[j][i];
        const Real  h      = L/2;
        const Real  dx     = delta.x;
        const Real  Peff   = (h*P_i + (dx-h) * P_ip)/dx;
        const Real  d2P    = P_psi - (Peff+Peff) + P_ip;
        const Real  fac    = dx/h;
        L2[j][im].x        = (P_i+P_i) - P_ip + fac*fac*d2P;
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
    
    const unit_t i     = J->upper;
    const unit_t im    = i-1;
    const unit_t ip    = i+1;
    const Real   Xi    = X[i];
    const Real   phi   = K->value - Xi; assert(phi>=0);
    const Real   psi   = Xi - J->value; assert(psi>=0);
    const Real   P_phi = K->pressure;
    const Real   P_psi = J->pressure;
    const Real   L     = phi+psi;
    const Real   scale = two_delta.x/L;
    L1[j][im].x = P_psi * scale;
    E1[j][ip].x = P_phi * scale;
    
    const Real P_i  = P[j][i];
    const Real h    = L/2;
    Real       Peff = P_i;
    if( phi > psi )
    {
        Peff = (h * P_i + (phi - h) * P_phi) / phi;
    }
    else
    {
        if(psi > phi )
        {
            Peff = (h*P_i + (psi-h) * P_psi) / psi;
        }
    }
    const Real fac = delta.x/h;
    const Real d2P = P_psi - (Peff+Peff) + P_phi;
    L2[j][im].x = 
    E2[j][ip].x = P_i + 0.5 * fac * fac * d2P;
    
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
        const Real jm     = j-1;
        const Real jp     = j+1;
        const Real P_jm   = P[jm][i];
        const Real Yj     = Y[j];
        const Real phi    = J->value - Yj; assert(phi>=0);
        const Real P_phi  = J->pressure;
        const Real dy     = delta.y;
        const Real L      = dy+phi;
        const Real sigma  = (P_phi - P_jm) / L;
        
        E1[jp][i].y       = P_jm + two_delta.y * sigma;
        
        const Real P_j  = P[j][i];
        const Real h    = L/2;
        const Real Peff = (h*P_j + (dy-h)*P_jm) / dy;
        const Real d2P  = P_jm - (Peff+Peff) + P_phi;
        const Real fac  = dy/h;
        E2[jp][i].y     = (P_j+P_j) - P_jm + fac*fac * d2P;
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
        const unit_t jp    = j+1;
        const unit_t jm    = j-1;
        const Real   P_jp  = P[jp][i];
        const Real   Yj    = Y[j];
        const Real   psi   = Yj - K->value; assert(psi>=0);
        const Real   P_psi = K->pressure;
        const Real   dy    = delta.y;
        const Real   L     = dy+psi;
        const Real   sigma = (P_jp - P_psi) /L;
        L1[jm][i].y = P_jp - two_delta.y * sigma;
        
        const Real   P_j  = P[j][i];
        const Real   h    = L/2;
        const Real   Peff = (h*P_j + (dy-h) * P_jp)/dy;
        const Real   d2P  = P_psi - (Peff+Peff) + P_jp;
        const Real   fac  = dy/h;
        L2[jm][i].y       = (P_j+P_j) - P_jp + fac*fac*d2P;
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
    
    const unit_t j  = J->upper;
    const unit_t jm = j-1;
    const unit_t jp = j+1;
    
    const Real Yj    = Y[j];
    const Real phi   = K->value - Yj; assert(phi>=0);
    const Real P_phi = K->pressure;
    const Real psi   = Yj - J->value; assert(psi>=0);
    const Real P_psi = J->pressure;
    const Real L     = phi+psi;
    
    const Real scale = two_delta.y/L;
    L1[jm][i].y = P_psi * scale;
    E1[jp][i].y = P_phi * scale;
    
    const Real P_j   = P[j][i];
    const Real h     = L/2;
    
    
    Real       Peff = P_j;
    if( phi > psi )
    {
        Peff = (h * P_j + (phi - h) * P_phi) / phi;
    }
    else
    {
        if(psi > phi )
        {
            Peff = (h*P_j + (psi-h) * P_psi) / psi;
        }
    }
    const Real fac = delta.y/h;
    const Real d2P = P_psi - (Peff+Peff) + P_phi;

    
    L2[jm][i].y = 
    E2[jp][i].y = P_j + 0.5 * fac * fac *d2P;
    
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
    pressurize_horz();
    pressurize_vert();
}



