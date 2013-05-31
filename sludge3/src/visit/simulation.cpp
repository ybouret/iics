#include "simulation.hpp"

Simulation:: ~Simulation() throw()
{
    
    
}


Simulation:: Simulation(const mpi   &MPI,
                        const Coord  N,
                        const Vertex Q) :
Workspace(MPI,N,Q),
VisIt::Simulation(MPI)
{
    
}


void Simulation:: get_meta_data(visit_handle &md) const
{
    //! append the mesh
    {
        VisIt_SimulationMetaData_addMesh(md, mesh_meta_data(mesh, "grid", par_size));
    }
    
}


visit_handle Simulation:: get_mesh( int domain, const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if( name == "grid" )
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
            
#if 0
            int minRealIndex[3] = { 1,         2 ,        0 };
            int maxRealIndex[3] = { int(X.width)-2, int(Y.width)-3, 0 };
            
            if(  par_rank < par_size - 1)
                maxRealIndex[1]++;
            
            VisIt_RectilinearMesh_setRealIndices(h, minRealIndex, maxRealIndex);
#endif
        }
    }
    return h;
    
}
