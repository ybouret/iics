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
    VertexArray   &E1;  //!< pressure when entering a bubble along x or y, first order
    VertexArray   &L1;  //!< pressure when leaving  a bubble along x or y, first order
    VertexArray   &E2;  //!< pressure when entering a bubble along x or y, second order
    VertexArray   &L2;  //!< pressure when entering a bubble along x or y, second order

    Array         &DeltaP; //!< store laplacian pressure to check...
    bool           right_wall; //!< default: false
    Real           P_user;     //!< default: 0.5
    
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
    /**
     \return 1 if success, 0 otherwise
     */
    int update_pressure(const mpi &MPI,
                        ColorType c,
                        const Real ftol);
    
    
    void compute_pressure( const mpi &MPI, const Real ftol );
    
    
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
