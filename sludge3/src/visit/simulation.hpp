#ifndef SLUDGE_VISIT_SIMULATION_INCLUDED
#define SLUDGE_VISIT_SIMULATION_INCLUDED 1

#include "../mpi/workspace.hpp"
#include "yocto/spade/visit.hpp"
#include "yocto/visit/interface.hpp"


class Simulation :
public Workspace,
public VisIt::Simulation,
public VisItIO
{
public:
    explicit Simulation(const mpi   &ref,
                        const Coord  N,
                        const Vertex Q
                        );
    virtual ~Simulation() throw();
    
    virtual void         get_meta_data(visit_handle &md ) const;
    virtual visit_handle get_mesh( int domain, const string &name ) const;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};

#endif

