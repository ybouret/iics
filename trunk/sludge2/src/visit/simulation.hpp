#ifndef SIM_INCLUDED
#define SIM_INCLUDED 1

#include "../mpi/cell.hpp"
#include "yocto/visit/interface.hpp"
#include "yocto/spade/visit.hpp"

class Simulation :
public VisIt::Simulation,
public Cell,
public VisItIO
{
public:
    explicit Simulation( const mpi &MPI );
    virtual ~Simulation() throw();
    
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};

#endif
