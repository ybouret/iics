#include "../cell.hpp"



void Cell:: compute_gradP()
{
    
    //--------------------------------------------------------------------------
    // initialize: ok for walls and inside bubbles
    //--------------------------------------------------------------------------
    gradP.ldz();
    
    const unit_t xlo   = lower.x + 1;
    const unit_t xhi   = upper.x;
    const unit_t xhim1 = xhi-1;
    
    for( unit_t j=lower.y;j<=upper.y;++j)
    {
        VertexArray1D & g_j = gradP[j];
        const Array1D & B_j = B[j];
        const Array1D & P_j = P[j];
        
        for( unit_t i=xlo;i<xhi;++i)
        {
            if(B_j[i]<=0)
            {
                const Real p_left  = P_left( j,i);
                const Real p_right = P_right(j,i);
                const Real p_lower = P_lower(j,i);
                const Real p_upper = P_upper(j,i);
                Vertex &g = g_j[i];
                g.x = inv_twodel.x * (p_right - p_left );
                g.y = inv_twodel.y * (p_upper - p_lower);
            }
        }
        
        //----------------------------------------------------------------------
        // process right boundary condition
        //----------------------------------------------------------------------
        if( right_wall )
        {
            //-- zero right gradient, already set
        }
        else
        {
            //-- right gradient estimation
            Vertex &g = g_j[xhi];
            g.y = inv_twodel.y * (P[j+1][xhi] -P[j-1][xhi]);
            if( B_j[xhim1] <= 0 )
            {
                //-- second order
                g.x = inv_twodel.x * ( 3*P_j[xhi] - 4* P_j[xhim1] + P_left(j,xhim1));
            }
            else
            {
                //-- first order only
                g.x = inv_delta.x * ( P_j[xhi] - P_left(j,xhi) );
            }
        }
        
    }
    
    
}
