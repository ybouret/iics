#include "simulation.hpp"

Simulation:: ~Simulation() throw()
{
}


Simulation:: Simulation(unit_t        Nx,
                        unit_t        Ny,
                        const Vertex &box,
                        const mpi    &mpi_ref) :
Cell(Nx,Ny,box,mpi_ref),
VisIt::Simulation(mpi_ref),
_visit()
{
}

void Simulation:: initialize()
{
    const double Ly = sim_region.length.y;
    const double Lx = sim_region.length.x;
    U.ldz();
    P.ldz();
    bubbles.empty();
    for( unit_t j=lower.y;j<=upper.y;++j)
    {
        for( unit_t i=lower.x;i<=upper.x;++i)
        {
            P[j][i]   = Y[j] / Ly;
            P[j][i]   = 0.1 * (0.5 - Alea());
            U[j][i].y = 0.04 + 0.08 * cos( numeric<double>::pi * Y[j] / Ly );
            U[j][i].y = 0.15;
        }
    }
    if( master )
    {
        Bubble *b = bubbles.create();
        b->map_peanut( Vertex(sim_box.x/2,0), 0.15 * Lx, 0.9 + 0.09 * Alea() );
        b->pressure = 1;
    }
    compute_fields();
}
