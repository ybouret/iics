#ifndef SLUDGE_MPI_WORKSPACE_INCLUDED
#define SLUDGE_MPI_WORKSPACE_INCLUDED 1

#include "./parameters.hpp"
#include "yocto/spade/mpi/workspace.hpp"
#include "yocto/comparator.hpp"

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
    Array         &W;   //!< Weight of the P[j][i] factor in the Poisson equation
    Array         &Bulk;       //!< number of bulk around another bulk
    Array         &DeltaP;     //!< store laplacian pressure to check...
    VertexArray   &V;          //!< velocities on the grid
    bool           right_wall; //!< default: false
    Real           P_user;     //!< default: 0.5
    const unit_t   jcoll_lo;   //!< lower index for marker collection
    const unit_t   jcoll_up;   //!< upper index for marker collection
    
    //! regularize bubbles, check boundaries, broadcast is_valid
    void validate_bubbles(const mpi &MPI);
    void broadcast_bubbles(const mpi &MPI);
    
    //! perform the segmentation
    /**
     fill the B field, compute the Bulk field
     and collect all markers
     */
    void segment();
    
    void save_markers() const;
    
    
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
                        ColorType c);
    
    
    //! Red-Black Gauss-Seidel until ftol
    /**
     then the gradient is computed and
     velocities are deduced,
     in the bulk and for each marker
     */
    void compute_pressure( const mpi &MPI );
    
    
    //! to debug
    void compute_laplacian();
    
    
       
    //! validate/compute pressure/velocities
    /**
     \return false is invalid bubbles
     */
    bool initialize( const mpi &MPI );
    
    
    //! move tracers associated to markers
    /**
     \return false if invalid bubbles are created
     */
    bool evolution( const mpi &MPI, Real dt );

    
    Vertex gradP_to_V( const Vertex &g ) const;
    
    void save_markers( const mpi &MPI ) const;
    
   
    struct LocalPressure
    {
        Vertex r; //!< position/delta position
        Real   P; //!< pressure in field
        Real   d; //!< distance to tracer
        Vertex g; //!< gradient
        
        static inline
        int CompareByVertex( const LocalPressure &lhs, const LocalPressure &rhs) throw()
        {
            return Vertex::lexicographical_compare(lhs.r, rhs.r);
        }
        
        static inline
        int CompareByDecreasingDistance( const LocalPressure &lhs, const LocalPressure &rhs) throw()
        {
            return __compare<Real>(rhs.d,lhs.d);
        }
        
        static inline
        int CompareByIncreasingDistance( const LocalPressure &lhs, const LocalPressure &rhs) throw()
        {
            return __compare<Real>(lhs.d,rhs.d);
        }
        
        
        
    };
    
    void collect_pressure( const Junction *J, LocalPressure lp[], size_t &n) const;
    
    
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
    
    void compute_velocities();
    
};

#endif
