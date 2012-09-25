#include "../simulation.hpp"

void Simulation:: get_meta_data(visit_handle &md) const
{

    //! append an UI command
    {
        visit_handle cmd = VISIT_INVALID_HANDLE;
        if(VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY)
        {
            VisIt_CommandMetaData_setName(cmd, "raz" );
            VisIt_SimulationMetaData_addGenericCommand(md, cmd);
        }
    }
    
    //! append the mesh
    {
        visit_handle mmd = mesh_meta_data(mesh, MeshName, par_size);
        VisIt_SimulationMetaData_addMesh(md, mmd);
    }
    
    //! append B on mesh
    {
        visit_handle vmd = variable_meta_data<Real>("B", MeshName);
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
    //! append P on mesh
    {
        visit_handle vmd = variable_meta_data<Real>("P", MeshName);
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
}
