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
    bool           is_valid;   //!< if bubbles are ok
    const Array1D &X;
    const Array1D &Y;
    Array         &P;
    Array         &B;
    VertexArray   &gradP;
    VertexArray   &Enter; //!< pressure when entering a bubble along x or y
    VertexArray  &Leave; //!< pressure when leaving  a bubble along x or y
    
    //! regularize bubbles, check boundaries, broadcast is_valid
    void validate_bubbles(const mpi &MPI);
    void broadcast_bubbles(const mpi &MPI);
    
    
    
    //! perform the segmentation in B field
    void segment();
    
    //! set pressure inside bubbles 
    void pressurize_bubbles();
    void compute_gradient(const mpi &MPI);
    
    //! set pressure to zero then pressurize
    void reset_pressure();
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Workspace);
    vector<Real> bpres; //!< bubble pressures
};

#endif