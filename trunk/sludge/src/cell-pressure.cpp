#include "cell.hpp"
#include "yocto/swamp/mpi.hpp"

void Cell:: compute_pressure()
{
    
    //==========================================================================
    //
    // we start from the bubble, with the B field initialized
    //
    //==========================================================================
    MPI.Printf0( stderr, "\t---> computing pressure\n");

    
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
            gradP[m->coord].ldz();
        }
    }
    
    for( unit_t j=upper.y;j>=lower.y;--j)
    {
        P[j][upper.x] = 0.5;
    }
    //--------------------------------------------------------------------------
    //
    // Red Black Gauss Seidel
    //
    //--------------------------------------------------------------------------
    static const unit_t shift_tab[2][2] =
    {
        {0,1},
        {1,0}
    };
    
    const Real ftol = 1e-5;
    
    for( size_t iter=1;;++iter)
    {
        //----------------------------------------------------------------------
        //
        // synchronize pressure fields
        //
        //----------------------------------------------------------------------
        _mpi::synchronize_one(P, MPI, *this, requests);
        
        //----------------------------------------------------------------------
        //
        // apply Red/Black Gauss Seidel
        //
        //----------------------------------------------------------------------
        const unit_t xmax = upper.x - 1;
        const unit_t xmin = lower.x + 1;
        
        int cvg = 1;
        for(int red=0; red <=1; ++red )
        {
            
            for(unit_t j=upper.y;j>=lower.y;--j)
            {
                Array1D       &P_j     = P[j];
                const Array1D &B_j     = B[j];
                VertexArray1D &gradP_j = gradP[j];
                const int      r       = (j&1) ? 1 : 0;
                const unit_t   shift   = shift_tab[red][r];
                
                //--------------------------------------------------------------
                //-- left gradient
                //--------------------------------------------------------------
                gradP_j[lower.x].ldz();
                
                //--------------------------------------------------------------
                // right gradient
                //--------------------------------------------------------------
                // if(0==shift)
                if( in_walls )
                {
                    gradP_j[upper.x].ldz();
                }
                else
                {
                    gradP_j[upper.x].y = 0;
                    gradP_j[upper.x].x = inv_two_dX * ( 3*P_j[upper.x] + P_j[upper.x-2] - 4*P_j[upper.x-1]);
                }
                
                for( unit_t i=xmax-shift;i>=xmin; i -= 2 )
                {
                    Vertex &gradP_ji   = gradP_j[i];
                    const Real P_j_im  = P_j[i-1];
                    const Real P_j_ip  = P_j[i+1];
                    const Real P_jm_i  = P[j-1][i];
                    const Real P_jp_i  = P[j+1][i];
                    gradP_ji.x         = (P_j_ip - P_j_im) * inv_two_dX;
                    gradP_ji.y         = (P_jp_i - P_jm_i) * inv_two_dY;
                    if( B_j[i] <= 0 )
                    {
                        //-- in the laponite
                        Real   &P_ji       = P_j[i];
                        const Real P0      = P_ji;
                        const Real twoP0   = P0+P0;
                        const Real LPx     = (P_j_im - twoP0 + P_j_ip ) * inv_dX2;
                        const Real LPy     = (P_jm_i - twoP0 + P_jp_i ) * inv_dY2;
                        const Real residue = LPx+LPy;
                        const Real delta_P = residue * stencil_w;
                        P_ji -= delta_P;
                        if( Fabs(delta_P) > ftol * Fabs(P_ji) )
                            cvg = 0;
                    }
                }
                P_j[lower.x] = (4 * P_j[lower.x+1] - P_j[lower.x+2])/3;
                if(in_walls)
                {
                    P_j[upper.x] = (4 * P_j[upper.x-1] - P_j[upper.x-2])/3;
                }
                else
                {
                    gradP_j[upper.x].y =  (P[j+1][upper.x] - P[j-1][upper.x]) * inv_two_dX;
                }
            }
        }
        
        //----------------------------------------------------------------------
        // use Allreduce to find out if every one has converged
        //----------------------------------------------------------------------
        int converged = 0;
        MPI.Allreduce(&cvg, &converged, 1, MPI_INT, MPI_SUM,MPI_COMM_WORLD);
        //MPI.Printf0(stderr, "@%6lu : #converged=%d\n", iter, converged);
        if( sim_size == converged)
        {
            MPI.Printf0( stderr, "\t---> converged\n");
            break;
        }
    }
    
    
    
}