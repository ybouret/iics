#include "simulation.hpp"
#include "yocto/swamp/mpi.hpp"


Simulation:: Simulation(unit_t     Nx, 
                        unit_t     Ny,
                        Real       Lx,
                        Real       Ly,
                        const mpi &ref) :
Cell(Nx,Ny,Lx,Ly,ref),
VisIt::Simulation(ref),
requests( num_requests())
{
    prepare_ghosts();
}

Simulation:: ~Simulation() throw()
{
    
}

void Simulation:: init_exchange()
{
    _mpi::init_exchange(MPI, *this, requests);
}

void Simulation:: wait_exchange()
{
    _mpi::wait_exchange(MPI, *this, requests);
}
