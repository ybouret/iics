#include "simulation.hpp"

visit_handle Simulation:: get_variable( int domain, const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    
    if( name == "P" )
    {
        const int nComponents= 1;
        const int nTuples    = P.items;
        //MPI.Printf0( stderr, "Sending P: %dx%d\n", nComponents, nTuples);
        assert(P.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, P.entry);
        }
        return h;
    }
    
    if( name == "B" )
    {
        const int nComponents= 1;
        const int nTuples    = B.items;
        //MPI.Printf0( stderr, "Sending B: %dx%d\n", nComponents, nTuples);
        assert(B.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, B.entry);
        }
        return h;
    }
    
    
    if( name == "U" )
    {
        const int nComponents= 2;
        const int nTuples    = U.items;
        //MPI.Printf0( stderr, "Sending U: %dx%d\n", nComponents, nTuples);
        assert(U.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real*)(U.entry));
        }
        return h;
    }
    
    if( name == "gradP" )
    {
        const int nComponents= 2;
        const int nTuples    = gradP.items;
        //MPI.Printf0( stderr, "Sending U: %dx%d\n", nComponents, nTuples);
        assert(gradP.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real*)(gradP.entry));
        }
        return h;
    }

    
    
    return h;
}


