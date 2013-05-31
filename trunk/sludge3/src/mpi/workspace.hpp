#ifndef SLUDGE_MPI_WORKSPACE_INCLUDED
#define SLUDGE_MPI_WORKSPACE_INCLUDED 1

#include "./parameters.hpp"
#include "yocto/spade/mpi/workspace.hpp"


typedef mpi_workspace<Layout,rmesh,Real> WorkspaceType;


class Workspace : public Parameters, public WorkspaceType
{
public:
    virtual ~Workspace() throw();
    explicit Workspace(const mpi    &MPI,
                       const Coord   N,
                       const Vertex  Q);
    
    Bubbles        bubbles;
    Junctions      junctions;
    const Array1D &X;
    const Array1D &Y;
    Array         &P;
    Array         &B;
    VertexArray   &gradP;
    
    void broadcast_bubbles(const mpi &MPI);
    void segment();
    
    //! set pressure inside bubbles and compute boundary pressure fields
    void pressurize();
    
    //! set pressure to zero then pressurize
    void reset_pressure();
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Workspace);
    vector<Real> bpres; //!< bubble pressures
};

#endif
