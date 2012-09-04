#include "cell.hpp"

void Cell:: compute_velocities()
{
    MPI.Printf0(stderr, "\t---> compute velocities\n");
    const unit_t xmax = upper.x;
    const unit_t xmin = lower.x + 1;
    
    //==========================================================================
    //
    // On the grid
    //
    //==========================================================================
    if( bubbles_velocities )
    {
        MPI.Printf0(stderr, "\t\t [[ computed ]]\n");
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
                    U_j[i] = velocity_from(gradP_j[i]);
                }
                else
                {
                    //--------------------------------------------------------------
                    // in a bubble
                    //--------------------------------------------------------------
                    U_j[i].ldz(); // idem if gradP_j[i] == 0.
                }
            }
        }
        
#if 0
        for( unit_t j=upper.y;j>=lower.y;--j)
        {
            for( unit_t i=upper.x;i>=lower.x;--i)
            {
                U[j][i].x = Real(i)/width.x;
                U[j][i].y = Real(j)/width.y;
            }
        }
#endif
        
    }
    else
    {
         MPI.Printf0(stderr, "\t\t [[ NOT computed ]]\n");
    }
    //==========================================================================
    //
    // On the bubbles
    //
    //==========================================================================
    
    for( Bubble *bubble=bubbles.first(); bubble; bubble=bubble->next)
    {
        for( Spot *spot = bubble->spots.head;spot;spot=spot->next)
        {
            const unit_t i0 = spot->gLower.x;
            const unit_t j0 = spot->gLower.y;
            const unit_t i1 = spot->gUpper.x;
            const unit_t j1 = spot->gUpper.y;
            const Real   x   = spot->bary.x;
            const Real   y   = spot->bary.y;
            const Real   umx = 1-x;
            const Real   umy = 1-y;
            
            Vertex      v;
            
            const Real  w00 = umx*umy;
            v += w00 * U[j0][i0];
            
            const Real  w10 = x*umy;
            v += w10 * U[j0][i1];
            
            const Real  w11 = x*y;
            v += w11 * U[j1][i1];
            
            const Real  w01 = umx * y;
            v += w01 * U[j1][i0];
            
            spot->has_U = true;
            //spot->U     = velocity_from(v);
            spot->U = v;
            
#if 0
            Vertex      u(0,0);
            Real  weight    = 0;
            
            if( B[j0][i0] <= 0 )
            {
                const Real  w00 = umx*umy;
                u += w00 * U[j0][i0];
                weight += w00;
            }
            
            if( B[j0][i1] <= 0)
            {
                const Real  w10 = x*umy;
                u += w10 * U[j0][i1];
                weight += w10;
            }
            
            if( B[j1][i1] <= 0)
            {
                const Real  w11 = x*y;
                u += w11 * U[j1][i1];
                weight += w11;
            }
            
            if( B[j1][i0] <= 0 )
            {
                const Real  w01 = umx * y;
                u += w01 * U[j1][i0];
                weight += w01;
            }
            
            if( weight > 0 )
            {
                spot->U     = (1/weight) * u;
                spot->has_U = true;
            }
            else
            {
                spot->has_U = false;
                spot->U.ldz();
            }
#endif
        }
    }
    
}