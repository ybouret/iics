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
    
    //==========================================================================
    //
    // synchronize gradP
    //
    //==========================================================================
    sync1(MPI, gradP);
    
}