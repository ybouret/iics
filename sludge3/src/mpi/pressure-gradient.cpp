#include "workspace.hpp"



void Workspace:: compute_gradP(const mpi &MPI)
{
    
    gradP.ldz();
    for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
    {
        for(unit_t i=bulk_imin;i<=bulk_imax;++i)
        {
            const Real which = B[j][i];
            if(which<0)
            {
                Vertex &g = gradP[j][i];
                g.x = (Enter[j][i+1].x - Leave[j][i-1].x) * order1fac.x;
                g.y = (Enter[j+1][i].y - Leave[j-1][i].y) * order1fac.y;
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