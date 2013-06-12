#include "simulation.hpp"

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
    
    if( name == "E1" )
    {
        const int nComponents= 2;
        const int nTuples    = E1.items;
        assert(E1.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real *)(E1.entry));
        }
    }
    
    if( name == "L1" )
    {
        const int nComponents= 2;
        const int nTuples    = L1.items;
        assert(L1.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real *)(L1.entry));
        }
    }
    
    if( name == "E2" )
    {
        const int nComponents= 2;
        const int nTuples    = E1.items;
        assert(E1.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real *)(E2.entry));
        }
    }
    
    if( name == "L2" )
    {
        const int nComponents= 2;
        const int nTuples    = L1.items;
        assert(L1.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real *)(L2.entry));
        }
    }
    
    
    if( name == "DeltaP" )
    {
        const int nComponents= 1;
        const int nTuples    = DeltaP.items;
        assert(DeltaP.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, DeltaP.entry);
        }
    }
    
    if( name == "W" )
    {
        const int nComponents= 1;
        const int nTuples    = W.items;
        assert(W.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, W.entry);
        }
    }
    
    if( name == "Bulk" )
    {
        const int nComponents= 1;
        const int nTuples    = Bulk.items;
        assert(Bulk.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, Bulk.entry);
        }
    }
    
    
    
    if( name == "V" )
    {
        const int nComponents= 2;
        const int nTuples    = V.items;
        assert(V.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real*)V.entry);
        }
    }
    
    
    return h;
    
}
