#include "../simulation.hpp"

visit_handle Simulation:: get_mesh( int domain, const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if( name == MeshName )
    {
        if( VisIt_RectilinearMesh_alloc(&h) == VISIT_OKAY )
        {
            assert( h != VISIT_INVALID_HANDLE );
            
            visit_handle   hx = VISIT_INVALID_HANDLE;
            visit_handle   hy = VISIT_INVALID_HANDLE;
            VisIt_VariableData_alloc(&hx);
            VisIt_VariableData_setDataD(hx, VISIT_OWNER_SIM, 1, X.items, X.entry );
            VisIt_VariableData_alloc(&hy);
            VisIt_VariableData_setDataD(hy, VISIT_OWNER_SIM, 1, Y.items, Y.entry );
            
            VisIt_RectilinearMesh_setCoordsXY(h, hx, hy);
            
            int minRealIndex[3] = { 0,         2 ,        0 };
            int maxRealIndex[3] = { X.width-1, Y.width-2, 0 };
                        
            VisIt_RectilinearMesh_setRealIndices(h, minRealIndex, maxRealIndex);
        }
    }
    return h;
}

