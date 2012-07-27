#include "cell.hpp"

void Cell:: compute_velocities()
{
    MPI.Printf0(stderr, "\t---> compute velocities\n");
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
    
    for( Bubble *bubble=bubbles.first(); bubble; bubble=bubble->next)
    {
        for( Spot *spot = bubble->spots.head;spot;spot=spot->next)
        {
            //Tracer      *p  = spot->handle;
            const unit_t i0 = spot->gLower.x;
            const unit_t j0 = spot->gLower.y;
            const unit_t i1 = spot->gUpper.x;
            const unit_t j1 = spot->gUpper.y;
            
            //------------------------------------------------------------------
            // interpolate the gradient
            //------------------------------------------------------------------
            Vertex       g(0,0);
            const Real  x   = spot->bw.x;
            const Real  y   = spot->bw.y;
            //fprintf(stderr, "wx=%.4f,wy=%.4f\n", x,y);
            const Real  umx = 1-x;
            const Real  umy = 1-y;
            const Real  w00 = umx*umy;
            const Real  w10 = x*umy;
            const Real  w11 = x*y;
            const Real  w01 = umx * y;
            
            
            g += w00 * gradP[j0][i0];
            g += w10 * gradP[j0][i1];
            g += w11 * gradP[j1][i1];
            g += w01 * gradP[j1][i0];

            
                        
            spot->U = velocity_from(g);
        }
    }
    
}