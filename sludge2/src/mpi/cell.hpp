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
                  const Vertex      &L,
                  const FieldsSetup &F);
    
    virtual ~Cell() throw();
    const Array1D &X;
    const Array1D &Y;
    
    Segmenter  segmenter;
    Bubbles    bubbles;
    const Real ymin;
    const Real ymax;
    
    //! broadcast and find spots
    void dispatch( const mpi &MPI );
    
   
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};

#endif

