#include "../cell.hpp"

void Cell:: compute_gradP()
{
    
    //--------------------------------------------------------------------------
    // initialize
    //--------------------------------------------------------------------------
    gradP.ldz();
    
    
    //--------------------------------------------------------------------------
    // compute where it is possible
    //--------------------------------------------------------------------------
    for( unit_t j=lower.y;j<=upper.y;++j)
    {
        VertexArray1D &g_j = gradP[j];
        for(unit_t i=lower.x+1;i<upper.x;++i)
        {
            if( B[j-1][i] <=0 && B[j+1][i] <=0 && B[j][i-1]<=0 && B[j][i+1] <=0 )
            {
                g_j[i].x = inv_twodel.x*(P[j][i+1]-P[j][i-1]);
                g_j[i].y = inv_twodel.y*(P[j+1][i]-P[j-1][i]);
            }
        }
    }
    
    
}
