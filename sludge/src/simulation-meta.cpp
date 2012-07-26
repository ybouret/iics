#include "simulation.hpp"

const char Simulation::MeshName[] = "mesh";

void Simulation:: get_meta_data( visit_handle &md ) const
{
    
    //! append an UI command
    {
        visit_handle cmd = VISIT_INVALID_HANDLE;
        if(VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY)
        {
            VisIt_CommandMetaData_setName(cmd, "save1" );
            VisIt_SimulationMetaData_addGenericCommand(md, cmd);
        }
    }
    
    
    
    //! append an UI command
    {
        visit_handle cmd = VISIT_INVALID_HANDLE;
        if(VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY)
        {
            VisIt_CommandMetaData_setName(cmd, "raz" );
            VisIt_SimulationMetaData_addGenericCommand(md, cmd);
        }
    }
    
    {
        visit_handle cmd = VISIT_INVALID_HANDLE;
        if(VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY)
        {
            VisIt_CommandMetaData_setName(cmd, "inter" );
            VisIt_SimulationMetaData_addGenericCommand(md, cmd);
        }
    }
    
    //! append the mesh
    {
        visit_handle mmd = _visit::mesh_meta_data(mesh, MeshName, par_size);
        VisIt_SimulationMetaData_addMesh(md, mmd);
    }
    
    if( !parallel )    //! append the bubbles ?
    {
        const Bubble *p = bubbles.first();
        while( p )
        {
            visit_handle cmd = VISIT_INVALID_HANDLE;
            if( VisIt_CurveMetaData_alloc(&cmd) == VISIT_OKAY )
            {
                const string bubble_name = "bubble";
                VisIt_CurveMetaData_setName(cmd, bubble_name.c_str());
                
                VisIt_SimulationMetaData_addCurve(md, cmd);
            }
            p=p->next;
            break;
        }
    }
    
    //! append P on mesh
    {
        visit_handle vmd = variable_meta_data<Real>("P", MeshName);
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
    //! append B on mesh
    {
        visit_handle vmd = variable_meta_data<Real>("B", MeshName);
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
    //! append U on mesh
    {
        visit_handle vmd = variable_meta_data<Vertex>("U", MeshName);
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
    //! append gradP on mesh
    {
        visit_handle vmd = variable_meta_data<Vertex>("gradP", MeshName);
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
    
}
