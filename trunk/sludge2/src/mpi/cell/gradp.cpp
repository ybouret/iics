#include "../cell.hpp"



void Cell:: compute_gradP()
{
    
    //--------------------------------------------------------------------------
    // initialize: ok for walls and inside bubbles
    //--------------------------------------------------------------------------
    gradP.ldz();
    
    const unit_t xlo = lower.x + 1;
    const unit_t xhi = upper.x;
    for( unit_t j=lower.y;j<=upper.y;++j)
    {
        VertexArray1D & g_j = gradP[j];
        const Array1D & B_j = B[j];
        
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
        
    }
    
    
}
