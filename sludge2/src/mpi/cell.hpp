#ifndef CELL_INCLUDED
#define CELL_INCLUDED

#include "../parameters.hpp"
#include "yocto/spade/mpi/workspace.hpp"

typedef mpi_workspace<Layout,rmesh,Real> Workspace;

//! Hele Shaw Cell
class Cell : public Parameters, public Workspace
{
public:
    explicit Cell(const mpi         &MPI,
                  const Coord       &N,
                  const Vertex      &L);
    
    virtual ~Cell() throw();
    const Array1D &X;
    const Array1D &Y;
    Array         &B;     //!< bubble markers
    Array         &P;     //!< pressure
    VertexArray   &gradP; //!< pressure gradient
    VertexArray   &U;     //!< velocity field
    Segmenter      segmenter;
    Bubbles        bubbles;
    const Real     ymin;
    const Real     ymax;
    
    //! broadcast, find spots and bubbles, sync B field
    void dispatch( const mpi &MPI );
    
    //! save gnuplot B field
    void save_B( const string &filename ) const;
    
    //! compute pressure from bubbles
    void compute_pressure( const mpi &MPI );
    
    //! logical serie of event
    /**
     - master update topology
     - dispatch
     - compute pressure
     */
    void legalize( const mpi &MPI );
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};

#endif

