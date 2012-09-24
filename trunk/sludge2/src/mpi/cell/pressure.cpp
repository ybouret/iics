#include "../cell.hpp"

void Cell:: init_pressure(const mpi &MPI )
{
    //! bubble
    segmenter.pressurize(P);
    
    //! boundary conditions
    
    //! sync pressure
    sync1(MPI,P);
    
}