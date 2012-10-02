#ifndef CELL_INCLUDED
#define CELL_INCLUDED

#include "../parameters.hpp"
#include "yocto/spade/mpi/workspace.hpp"
#include "yocto/spade/vtk/writer.hpp"

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
    Array         &B;         //!< bubble markers
    Array         &P;         //!< pressure
    VertexArray   &gradP;     //!< pressure gradient
    VertexArray   &U;         //!< velocity field
    VertexArray   &Penter;    //!< effective pressure in bubble, along X and along Y
    VertexArray   &Pleave;    //!< effective pressure when leaving bubble
    Segmenter      segmenter;
    Bubbles        bubbles;
    const Real     ymin;
    const Real     ymax;
    vtk_writer     vtk;
    
    //! broadcast, find spots and bubbles, sync B field
    void dispatch( const mpi &MPI );
    
    //! save gnuplot B field
    void save_B( const string &filename ) const;
    
    //! save gnuplot B field, for outline
    void save_outB( const string &filename ) const;
    
    //! compute P gradient from bubbles + boundary conditions
    void compute_gradP();
    
    //! compute pressure from bubbles
    /**
     assume bubbles are pressurized by segmenter !!
     ie that the fields Penter et Pleave
     are computed and sync.
     */
    void compute_pressure( const mpi &MPI );
    
    //! logical serie of event
    /**
     - master update topology
     - dispatch
     - location
     - segmentation
     - build effective pressure and fill bubble pressure
     - call compute pressure
     */
    void legalize( const mpi &MPI );
    
    Real P_left(  unit_t j, unit_t i) const throw();  //!< at j,i-1: assuming j,i in the bulk
    Real P_right( unit_t j, unit_t i) const throw();  //!< at j,i+1: assuming j,i in the bulk
    Real P_lower( unit_t j, unit_t i) const throw();  //!< at j-1,i: assuming j,i in the bulk
    Real P_upper( unit_t j, unit_t i) const throw();  //!< at j+1,i: assuming j,i in the bulk

    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};

#endif

