#ifndef SIMULATION_INCLUDED
#define SIMULATION_INCLUDED 1

#include "cell.hpp"
#include "yocto/swamp/visit.hpp"
#include "yocto/visit/interface.hpp"

class Simulation :
public Cell,
public VisIt::Simulation,
public _visit
{
public:
    Simulation(unit_t        Nx,
               unit_t        Ny,
               const Vertex &box,
               const mpi    &mpi_ref);
    
    virtual ~Simulation() throw();
    
    
    //! create initial bubble(s) and set fields
    /**
     no ghosts/bubble exchange
     */
    void initialize();
    
    virtual void         get_meta_data( visit_handle &md ) const;
    virtual visit_handle get_mesh( int domain, const string &name ) const;
    virtual visit_handle get_variable( int domain, const string &name ) const;
    virtual visit_handle get_curve( const string &name ) const;
    virtual void         perform( const string &cmd );
    
    virtual void step();
    
    static const char MeshName[];
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
    
};

#endif
