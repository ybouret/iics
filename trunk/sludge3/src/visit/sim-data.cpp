#include "simulation.hpp"

visit_handle Simulation:: get_variable( int domain, const string &name ) const
{
    
    if(name == "P")      return variable_data(P);
    if(name == "B")      return variable_data(B);
    if(name == "gradP")  return variable_data(gradP);
    if(name == "E1")     return variable_data(E1);
    if(name == "L1")     return variable_data(L1);
    if(name == "E2")     return variable_data(E2);
    if(name == "L2")     return variable_data(L2);
    if(name == "DeltaP") return variable_data(DeltaP);
    if(name == "W")      return variable_data(W);
    if(name == "Bulk")   return variable_data(Bulk);
    if(name == "V")      return variable_data(V);
    
    return VISIT_INVALID_HANDLE;
    
}
