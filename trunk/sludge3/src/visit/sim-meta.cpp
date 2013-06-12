#include "simulation.hpp"

void Simulation:: get_meta_data(visit_handle &md) const
{
    
    //! append an UI generic command
    add_generic_command("raz", md);
    
    //! append the grid
    VisIt_SimulationMetaData_addMesh(md, mesh_meta_data(mesh, "grid", par_size));
    
    
    //! append P on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Real>("P", "grid"));
    
    //! append V on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("V", "grid"));
    
    
    //! append B on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Real>("B", "grid"));
    
    //! append gradP on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("gradP", "grid"));
    
    //! append Enter on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("E1", "grid"));
    
    //! append Leave on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("L1", "grid"));
    
    //! append Enter on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("E2", "grid"));
    
    //! append Leave on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("L2", "grid"));
    
    //! append W on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Real>("W", "grid"));
    
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Real>("Bulk", "grid"));
    
    
    //! append DeltaP on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Real>("DeltaP", "grid"));
    
    //! append bubbles
    for( const Bubble *b = bubbles.head;b;b=b->next)
    {
        const string bn = vformat("bubble%u", unsigned(b->UID));
        visit_handle h = VISIT_INVALID_HANDLE;
        if(VisIt_CurveMetaData_alloc(&h) == VISIT_OKAY )
        {
            VisIt_CurveMetaData_setName(h, bn.c_str());
            //VisIt_CurveMetaData_setXLabel(cmd, "Angle");
            ////VisIt_CurveMetaData_setXUnits(cmd, "radians");
            //VisIt_CurveMetaData_setYLabel(cmd, "Amplitude");
            VisIt_SimulationMetaData_addCurve(md, h);
        }
    }
    
    //! append junctions
    if( !parallel )
    {
        visit_handle h = VISIT_INVALID_HANDLE;
        if(VisIt_CurveMetaData_alloc(&h) == VISIT_OKAY )
        {
            VisIt_CurveMetaData_setName(h, "junctions" );
            VisIt_SimulationMetaData_addCurve(md, h);
        }
    }
    
}
