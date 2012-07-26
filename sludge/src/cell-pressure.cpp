#include "cell.hpp"
#include "yocto/swamp/mpi.hpp"

void Cell:: compute_pressure()
{
    
    //==========================================================================
    //
    // we start from the bubble, with the B field initialized
    //
    //==========================================================================
    
    
    //--------------------------------------------------------------------------
    //
    // initialize pressure field
    //
    //--------------------------------------------------------------------------
    
#if 0
    //! debug only
    for( unit_t j=lower.y;j<=upper.y;++j)
    {
        for( unit_t i=lower.x;i<=upper.x;++i)
        {
            P[j][i] = (0.5-Alea());
        }
    }
#endif
    
    for( const Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next)
    {
        const Real pressure = bubble->pressure;
        for( const Marker *m = bubble->markers.head;m;m=m->next)
        {
            P[m->coord] = pressure;
        }
    }
    
    for( size_t iter=1; iter<=100;++iter)
    {
        //----------------------------------------------------------------------
        //
        // synchronize pressure fields
        //
        //----------------------------------------------------------------------
        _mpi::synchronize1(P, MPI, *this, requests);
        
        //----------------------------------------------------------------------
        //
        // apply Red/Black Gauss Seidel
        //
        //----------------------------------------------------------------------
        const unit_t xmax = upper.x - 1;
        const unit_t xmin = lower.x + 1;
        for(int r=0; r <=1; ++r )
        {
            for(unit_t j=upper.y;j>=lower.y;--j)
            {
                Array1D       &P_j = P[j];
                const Array1D &B_j = B[j];
                for( unit_t i=xmax-r;i>=xmin; i -= 2 )
                {
                    if( B_j[i] <= 0 )
                    {
                        //-- in the laponite
                        Real &P_ji         = P_j[i];
                        const Real P0      = P_ji;
                        const Real twoP0   = P0+P0;
                        const Real LPx     = (P_j[i-1]  - twoP0 + P_j[i+1] ) * inv_dX2;
                        const Real LPy     = (P[j-1][i] - twoP0 + P[j+1][i]) * inv_dY2;
                        const Real residue = LPx+LPy;
                        P_ji -= residue * stencil_w;
                    }
                }
            }
        }
    }
    
    
    
    
}