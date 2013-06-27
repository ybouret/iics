#include "workspace.hpp"



void Workspace:: compute_gradP(const mpi &MPI)
{
    
    gradP.ldz();
    for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
    {
        const unit_t jm = j-1;
        const unit_t jp = j+1;
        
        for(unit_t i=bulk_imin;i<=bulk_imax;++i)
        {
            const Real which = B[j][i];
            if(which<0)
            {
                Vertex &g = gradP[j][i];
                g.x = (E1[j][i+1].x - L1[j][i-1].x) * order1fac.x;
                g.y = (E1[jp][i].y  - L1[jm][i].y ) * order1fac.y;
            }
        }
    }
    
    if( !right_wall )
    {
        const unit_t i   = upper.x;
        const unit_t im  = i-1;
        const unit_t im2 = i-2;
        
        for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
        {
            const unit_t jm = j-1;
            const unit_t jp = j+1;
            Vertex &g = gradP[j][upper.x];
            g.y = (E1[jp][i].y  - L1[jm][i].y ) * order1fac.y;
            if(B[j][im]<0)
            {
                g.x = (4*P[j][im] - 3*P[j][i] - L1[j][im2].x) * order1fac.x;
            }
            else
            {
                g.x = (L1[j][im].x - P[j][i]) / delta.x;
            }
        }
    }
    
    //==========================================================================
    //
    // synchronize gradP
    //
    //==========================================================================
    sync1(MPI, gradP);
    
    
}