#include "cell.hpp"

void Cell:: compute_velocities()
{
    const unit_t xmax = upper.x - 1;
    const unit_t xmin = lower.x + 1;
    
    for( unit_t j=upper.y;j>=lower.y;--j)
    {
        const Array1D       &B_j     = B[j];
        const VertexArray1D &gradP_j = gradP[j];
        VertexArray1D       &U_j     = U[j];
        
        U_j[lower.x].ldz();
        
        for(unit_t i=xmax;i>=xmin;--i)
        {
            if(B_j[i]<=0)
            {
                //--------------------------------------------------------------
                // in the laponite
                //--------------------------------------------------------------
                U_j[i] = -gradP_j[i];
            }
            else
            {
                //--------------------------------------------------------------
                // in a bubble
                //--------------------------------------------------------------
                U_j[i].ldz();
            }
        }
    }
    
}