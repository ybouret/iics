#include "workspace.hpp"

#define BUBBLE_NONE    (0)
#define BUBBLE_BEFORE  (1)
#define BUBBLE_AFTER   (2)
#define BUBBLE_BOTH   ( BUBBLE_BEFORE | BUBBLE_AFTER )

void Workspace:: pressurize_horz()
{
    //==========================================================================
    //
    // using horizontal axis => x components of Leave/Enter
    //
    //==========================================================================
    for(unit_t j=outline.lower.y;j<=outline.upper.y;++j)
    {
        const Junction::List &JL = junctions.Horz(j);
        if(JL.size>0)
        {
            //==================================================================
            // There are some bubbles on the line
            //==================================================================
            
            //------------------------------------------------------------------
            //
            // in bulk
            //
            //------------------------------------------------------------------
            for(unit_t i=bulk_imin;i<=bulk_imax;++i)
            {
                if(B[j][i]<0)
                {
                    
                    //----------------------------------------------------------
                    // detect topology
                    //----------------------------------------------------------
                    const Real      Xi     = X[i];
                    int             status = BUBBLE_NONE;
                    const Junction *before = 0;
                    const Junction *after  = 0;
                    Real            phi    = 0;
                    Real            psi    = 0;
                    if(B[j][i-1]>=0)
                    {
                        status |= BUBBLE_BEFORE;
                        before  = JL.before(Xi);
                        if(!before)
                            throw exception("@y=%g: no junction before x=%g", Y[j], Xi);
                        psi = Xi - before->value;
                    }
                    
                    if(B[j][i+1]>=0)
                    {
                        status |= BUBBLE_AFTER;
                        after   = JL.after(Xi);
                        if(!after)
                            throw exception("@y=%g: no junction after x=%g", Y[j], Xi);
                        phi = after->value - Xi;
                    }
                    
                    //----------------------------------------------------------
                    // handle topology
                    //----------------------------------------------------------
                    switch(status)
                    {
                        case BUBBLE_BEFORE:
                        {
                            const Real Pp = P[j][i+1];
                            const Real Pb = before->pressure;
                            Leave[j][i-1].x = Pp - (2*delta.x) * (Pp - Pb)/(delta.x+psi);
                        }
                            break;
                            
                        case BUBBLE_AFTER:
                        {
                            const Real Pm   = P[j][i-1];
                            const Real Pb   = after->pressure;
                            Enter[j][i+1].x = Pm + (2*delta.x) * (Pb - Pm) / (delta.x + phi);
                        }
                            break;
                            
                        case BUBBLE_BOTH:
                        {
                            const Real fac = (2*delta.x)/(phi+psi);
                            Enter[j][i+1].x = after->pressure  * fac;
                            Leave[j][i-1].x = before->pressure * fac;
                            
                        }
                            break;
                            
                        default:
                            break;
                    }
                }
            }
            
            
            //------------------------------------------------------------------
            //
            // on Horizonal side(s)
            //
            //------------------------------------------------------------------
        }
        
    }
    
}

void Workspace:: pressurize_vert()
{
    //==========================================================================
    //
    // using Vertical axis => y components of Leave/Enter
    //
    //==========================================================================
    for(unit_t i=outline.lower.x; i <= outline.upper.x; ++i )
    {
        const Junction::List &JL = junctions.Vert(i);
        if(JL.size)
        {
            //==================================================================
            // There are some bubbles on the line
            //==================================================================
            
            //------------------------------------------------------------------
            //
            // in bulk
            //
            //------------------------------------------------------------------
            for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
            {
                if(B[j][i]<0)
                {
                    
                    //----------------------------------------------------------
                    // detect topology
                    //----------------------------------------------------------
                    const Real      Yj     = Y[j];
                    int             status = BUBBLE_NONE;
                    const Junction *before = 0;
                    const Junction *after  = 0;
                    Real            phi    = 0;
                    Real            psi    = 0;
                    if(B[j-1][i]>=0)
                    {
                        status |= BUBBLE_BEFORE;
                        before  = JL.before(Yj);
                        if(!before)
                            throw exception("@x=%g: no junction before y=%g", X[i], Yj);
                        psi = Yj - before->value;
                    }
                    
                    if(B[j+1][i]>=0)
                    {
                        status |= BUBBLE_AFTER;
                        after   = JL.after(Yj);
                        if(!after)
                            throw exception("@x=%g: no junction after y=%g", X[i], Yj);
                        phi = after->value - Yj;
                    }
                    
                    //----------------------------------------------------------
                    // handle topology
                    //----------------------------------------------------------
                    switch(status)
                    {
                        case BUBBLE_BEFORE:
                        {
                            const Real Pp = P[j+1][i];
                            const Real Pb = before->pressure;
                            Leave[j-1][i].y = Pp - (2*delta.y) * (Pp - Pb)/(delta.x+psi);
                        }
                            break;
                            
                        case BUBBLE_AFTER:
                        {
                            const Real Pm   = P[j-1][i];
                            const Real Pb   = after->pressure;
                            Enter[j+1][i].y = Pm + (2*delta.y) * (Pb - Pm) / (delta.x + phi);
                        }
                            break;
                            
                        case BUBBLE_BOTH:
                        {
                            const Real fac = (2*delta.y)/(phi+psi);
                            Enter[j+1][i].y = after->pressure  * fac;
                            Leave[j-1][i].y = before->pressure * fac;
                            
                        }
                            break;
                            
                        default:
                            break;
                    }
                    
                    
                }
            }
            
            //------------------------------------------------------------------
            //
            // on Vertical side(s)
            //
            //------------------------------------------------------------------
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



