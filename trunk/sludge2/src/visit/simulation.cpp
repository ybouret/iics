#include "simulation.hpp"


Simulation:: ~Simulation() throw()
{
}


Simulation:: Simulation( const mpi         &MPI,
                        const Coord       &N,
                        const Vertex      &L  ) :
VisIt::Simulation(MPI),
Cell(MPI,N,L),
VisItIO()
{

}


const char Simulation:: MeshName[] = "mesh";