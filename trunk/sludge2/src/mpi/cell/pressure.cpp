#include "../cell.hpp"

void Cell:: init_pressure(const mpi &MPI )
{
    //! pressure from bubbles
    segmenter.pressurize(P);
    
    //! boundary conditions
    for( unit_t j=lower.y;j<=upper.y;++j)
    {
        P[j][upper.x] = 0.5;
    }
    
    //! sync pressure
    sync1(MPI,P);
    
}