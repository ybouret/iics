#include "workspace.hpp"

void Workspace:: compute_velocities()
{
    for(unit_t j=outline.lower.y;j<=outline.upper.y;++j)
    {
        for(unit_t i=outline.lower.x;i<=outline.upper.x;++i)
        {
            V[j][i] = -gradP[j][i];
        }
        
    }
    
}
