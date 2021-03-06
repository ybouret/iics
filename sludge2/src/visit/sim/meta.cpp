#include "../simulation.hpp"

void Simulation:: get_meta_data(visit_handle &md) const
{
    
    add_generic_command( "raz", md);
    if( par_size == 1)
        add_generic_command( "geo", md);
    
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
    
    //! append Bulk on mesh
    {
        visit_handle vmd = variable_meta_data<Real>("Bulk", MeshName);
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
    //! append P on mesh
    {
        visit_handle vmd = variable_meta_data<Real>("P", MeshName);
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
    
    //! append Penter on mesh
    {
        visit_handle vmd = variable_meta_data<Vertex>("Penter", MeshName);
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
    //! append Pleave on mesh
    {
        visit_handle vmd = variable_meta_data<Vertex>("Pleave", MeshName);
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
    
    if( !parallel )    //! append the bubbles ?
    {
        const Bubble *p = bubbles.first();
        while(p)
        {
            visit_handle cmd = VISIT_INVALID_HANDLE;
            if( VisIt_CurveMetaData_alloc(&cmd) == VISIT_OKAY )
            {
                const string bubble_name = vformat("bubble%u",p->id);
                VisIt_CurveMetaData_setName(cmd, bubble_name.c_str());
                
                VisIt_SimulationMetaData_addCurve(md, cmd);
            }
            p=p->next;
            //break;
        }
        
        {
            visit_handle cmd = VISIT_INVALID_HANDLE;
            if( VisIt_CurveMetaData_alloc(&cmd) == VISIT_OKAY )
            {
                const string bubble_name = "junctions";
                VisIt_CurveMetaData_setName(cmd, bubble_name.c_str());
                
                VisIt_SimulationMetaData_addCurve(md, cmd);
            }
            
        }
    }
    
    
}
