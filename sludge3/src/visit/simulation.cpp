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
    
    //! append an UI generic command
    add_generic_command("raz1", md);
    
    //! append the grid
    VisIt_SimulationMetaData_addMesh(md, mesh_meta_data(mesh, "grid", par_size));
    
    
    //! append P on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Real>("P", "grid"));
    
    
    //! append B on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Real>("B", "grid"));

    //! append gradP on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("gradP", "grid"));

    
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


visit_handle Simulation:: get_variable( int domain, const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    
    if( name == "P" )
    {
        const int nComponents= 1;
        const int nTuples    = P.items;
        assert(P.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, P.entry);
        }
    }
    
    if( name == "B" )
    {
        const int nComponents= 1;
        const int nTuples    = B.items;
        assert(B.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, B.entry);
        }
    }
    
    if( name == "gradP" )
    {
        const int nComponents= 2;
        const int nTuples    = gradP.items;
        assert(gradP.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real *)(gradP.entry));
        }
    }
    
    return h;

}


void Simulation:: perform( const string &cmd, const array<string> &args)
{
    if( cmd == "raz1" )
    {
        init_one_bubble();
    }
}
