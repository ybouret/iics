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
    explicit Simulation(const mpi         &MPI,
                        const Coord       &N,
                        const Vertex      &L );
    virtual ~Simulation() throw();
    
    void initialize();
    
    //! register all meta data
    virtual void         get_meta_data( visit_handle &md ) const;
    
    //! get the mesh
    virtual visit_handle get_mesh( int domain, const string &name ) const;

    //! get data for mesh
    visit_handle get_variable( int domain, const string &name ) const;

    //! step
    virtual void step();
    
    static const char MeshName[];
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};

#endif
