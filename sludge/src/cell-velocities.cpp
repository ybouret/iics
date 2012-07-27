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
    
    for( Bubble *bubble=bubbles.first(); bubble; bubble=bubble->next)
    {
        for( Spot *spot = bubble->spots.head;spot;spot=spot->next)
        {
            Tracer      *p  = spot->handle;
            const unit_t i0 = p->gLower.x;
            const unit_t j0 = p->gLower.y;
            const unit_t i1 = p->gUpper.x;
            const unit_t j1 = p->gUpper.y;
            Vertex      u(0,0);
            const Real  x   = p->bw.x;
            const Real  y   = p->bw.y;
            const Real  umx = 1-x;
            const Real  umy = 1-y;
            Real        weight = 0;
            
            
            if(B[j0][i0]<=0)
            {
                const Real  w00 = umx*umy;
                weight += w00;
                u += w00 * U[j0][i0];
            }
            
            if(B[j0][i1]<=0)
            {
                const Real  w10 = x*umy;
                weight += w10;
                u += w10 * U[j0][i1];
            }
            
            if(B[j1][i1]<=0)
            {
                const Real  w11 = x*y;
                weight += w11;
                u += w11 * U[j1][i1];
            }
            
            if(B[j1][i0]<=0)
            {
                const Real  w01 = umx * y;
                weight += w01;
                u += w01 * U[j1][i0];
            }
            
            if( weight <= 0)
                throw exception("Invalid tracer for velocity!");
            
            
            spot->U = (1/weight) * u;
        }
    }
    
}