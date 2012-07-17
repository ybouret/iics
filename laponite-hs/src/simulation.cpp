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


void Simulation:: initialize()
{
    U.ldz();
    P.ldz();
    bubbles.none();
    for( unit_t j=lower.y;j<=upper.y;++j)
    {
        for( unit_t i=lower.x;i<=upper.x;++i)
        {
            P[j][i]   = Y[j] / Length.y;
            U[j][i].y = 0.04;
        }
    }
    if( master )
    {
        Bubble *b = bubbles.create();
        b->lambda = Length.x / (width.x*2);
        b->map_circle( V2D(Length.x/2,0), 0.2 * Length.y);
    }
}

void Simulation:: init_exchange()
{
    _mpi::init_exchange(MPI, *this, requests);
}

void Simulation:: wait_exchange()
{
    _mpi::wait_exchange(MPI, *this, requests);
}

void Simulation:: check_and_dispatch_bubbles()
{
    if( master )
    {
        master_update_topologies();
    }
    dispatch_bubbles(MPI);
}

