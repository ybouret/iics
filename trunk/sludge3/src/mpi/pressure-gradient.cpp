#include "workspace.hpp"



void Workspace:: compute_gradient(const mpi &MPI)
{

    gradP.ldz();
    for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
    {
        for(unit_t i=bulk_imin;i<=bulk_imax;++i)
        {
            const Real which = B[j][i];
            if(which>=0)
            {
                gradP[j][i].x = (P[j][i-1] - P[j][i+1]) * order1fac.x;
                gradP[j][i].y = (P[j+1][i] - P[j-1][i]) * order1fac.y;
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