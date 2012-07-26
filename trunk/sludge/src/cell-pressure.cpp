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
    for( const Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next)
    {
        const Real pressure = bubble->pressure;
        for( const Marker *m = bubble->markers.head;m;m=m->next)
        {
            P[m->coord] = pressure;
        }
    }
    
    //--------------------------------------------------------------------------
    //
    // Red Black Gauss Seidel
    //
    //--------------------------------------------------------------------------
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
                Array1D       &P_j     = P[j];
                const Array1D &B_j     = B[j];
                VertexArray1D &gradP_j = gradP[j];
                
                for( unit_t i=xmax-r;i>=xmin; i -= 2 )
                {
                    if( B_j[i] <= 0 )
                    {
                        //-- in the laponite
                        Vertex &gradP_ji   = gradP_j[i];
                        Real   &P_ji       = P_j[i];
                        const Real P0      = P_ji;
                        const Real twoP0   = P0+P0;
                        const Real P_j_im  = P_j[i-1];
                        const Real P_j_ip  = P_j[i+1];
                        const Real P_jm_i  = P[j-1][i];
                        const Real P_jp_i  = P[j+1][i];
                        gradP_ji.x         = (P_j_ip - P_j_im) * inv_two_dX;
                        gradP_ji.y         = (P_jp_i - P_jm_i) * inv_two_dY;
                        const Real LPx     = (P_j_im - twoP0 + P_j_ip ) * inv_dX2;
                        const Real LPy     = (P_jm_i - twoP0 + P_jp_i ) * inv_dY2;
                        const Real residue = LPx+LPy;
                        P_ji -= residue * stencil_w;
                    }
                }
            }
        }
    }
    
    // TODO: should be exchanged with U...
    _mpi::synchronize1(gradP, MPI, *this, requests);

    
    
}