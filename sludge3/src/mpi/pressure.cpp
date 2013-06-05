#include "workspace.hpp"


void Workspace:: update_pressure( const mpi &MPI, ColorType c )
{
    
    //--------------------------------------------------------------------------
    // set the boundary conditions
    //--------------------------------------------------------------------------
    for(unit_t j=outline.lower.x;j<=outline.upper.x;++j)
        P[j][upper.x] = 0.5;
    
    //--------------------------------------------------------------------------
    // prepare the pressure fields
    //--------------------------------------------------------------------------
    pressurize_bubbles();   //!< initialize Enter/Leave
    pressurize_contours();  //!< finalize   Enter/Leave according to bubbles
    
    //--------------------------------------------------------------------------
    // Evaluate gradient
    //--------------------------------------------------------------------------
    compute_gradP(MPI);
    
    //--------------------------------------------------------------------------
    // update according to color
    //--------------------------------------------------------------------------
    const unit_t shift[2] = { c, 1-c };
    for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
    {
        for( unit_t i=bulk_imin + shift[j&1]; i <= bulk_imax; i +=2 )
        {
            if(B[j][i]<0)
            {
                const Real P_center = P[j][i];
                const Real P_left   = Leave[j][i-1].x;
                const Real P_right  = Enter[j][i+1].x;
                const Real P_bottom = Leave[j-1][i].y;
                const Real P_top    = Enter[j+1][i].y;
                const Real mid      = -(P_center+P_center);
                DeltaP[j][i] = (P_left+mid+P_right) * order2fac.x + (P_bottom+mid+P_top) * order2fac.y;
            }
        }
    }
    
    
}
