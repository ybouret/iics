#include "cell.hpp"

void Cell:: compute_fields()
{
    dispatch_all();
    compute_pressure();
    compute_velocities();
    
}