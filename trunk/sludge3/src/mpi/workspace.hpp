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
    VertexArray   &Enter;  //!< pressure when entering a bubble along x or y
    VertexArray   &Leave;  //!< pressure when leaving  a bubble along x or y
    Array         &DeltaP; //!< store laplacian pressure to check...
    
    //! regularize bubbles, check boundaries, broadcast is_valid
    void validate_bubbles(const mpi &MPI);
    void broadcast_bubbles(const mpi &MPI);
    
    //! perform the segmentation in B field
    void segment();
    
    
    
    //! set pressure inside bubbles (once per step)
    void pressurize_bubbles();
    
    //! compute effective pressures using junctions topology
    void pressurize_contours();
    
    //! bubbles/contours must be pressurized
    void compute_gradP(const mpi &MPI);
    
    
    enum ColorType
    {
        Red=0,
        Black=1
    };
    
    //! assume bubbles are segmented
    void update_pressure( const mpi &MPI, ColorType c );
    
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Workspace);
    vector<Real> bpres; //!< bubble pressures
    void pressurize_horz();
    void pressurize_vert();
    
    void EnterX( const Junction *J, unit_t j);
    void LeaveX( const Junction *J, unit_t j);
    void AloneX( const Junction *J, const Junction *K, unit_t j);
    
    void EnterY( const Junction *J, unit_t i);
    void LeaveY( const Junction *J, unit_t i);
    void AloneY( const Junction *J, const Junction *K, unit_t i);
    
};

#endif
