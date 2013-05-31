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

