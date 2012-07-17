#ifndef SIMULATION_INCLUDED
#define SIMULATION_INCLUDED 1

#include "cell.hpp"
#include "yocto/visit/interface.hpp"

class Simulation : public Cell, public VisIt::Simulation 
{
public:
    Simulation(unit_t     Nx, 
               unit_t     Ny,
               Real       Lx,
               Real       Ly,
               const mpi &ref );
    virtual ~Simulation() throw();
    
    mpi::Requests requests;
    
    
    void init_exchange();
    void wait_exchange();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};

#endif
