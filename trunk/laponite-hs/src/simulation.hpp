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
    Simulation(unit_t     Nx, 
               unit_t     Ny,
               Real       Lx,
               Real       Ly,
               const mpi &ref );
    virtual ~Simulation() throw();
    
    //! requests for ghosts
    mpi::Requests requests;
    
    void check_and_dispatch_bubbles();
    
    //! start sending ghosts
    void init_exchange();
    
    //! finish sending ghosts
    void wait_exchange();
    
    
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
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Simulation);
};

#endif
