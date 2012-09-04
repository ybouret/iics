#ifndef CELL_INCLUDED
#define CELL_INCLUDED 1


#include "parameters.hpp"
#include "segmenter.hpp"
#include "rescaler.hpp"

//! Hele-Shaw Cell
class Cell : public Parameters, public WorkspaceBase
{
public:
    explicit Cell(unit_t        Nx, 
                  unit_t        Ny,
                  const Vertex &box,
                  const mpi    &mpi_ref);
    virtual ~Cell() throw();
    
    Array          &P;     //!< auxiliary, since it is synchronized alone
    Array          &B;     //!< bubble or not bubble, auxiliary
    VertexArray    &U;     //!< Ux,Uy
    VertexArray    &gradP; //!< pressure gradient
    const Array1D  &X;
    const Array1D  &Y;
    const Array1D  &dX;
    const Array1D  &dY;
    
    Bubbles            bubbles;
    Segmenter          segmenter;
    auto_ptr<Rescaler> rescaler;
    const Real         delta_X;
    const Real         delta_Y;
    const Real         inv_dX2;
    const Real         inv_dY2;
    const Real         inv_two_dX;
    const Real         inv_two_dY;
    const Real         stencil_w;       //!< inverse of P[i][j] factor in Laplacian
    Segment::List     *border_segments; //!< if parallel PBC
    const int          border_peer;     //!< with whom to complete segments
    const unit_t       border_j;        //!< indice
    const Real         border_y;        //!< to reconstruct the segments
    bool               in_walls;        //!< pressure boundary
    bool               bubbles_velocities;
    //==========================================================================
    // Bubbles Operations
    //==========================================================================
    //! broadcast bubbles, compute their properties and build segment them
    /**
     - bubbles.check_and_dispatch_all(MPI)
     - bubbles.check_geometries_within(sub_region.vmin.y, sub_region.vmax.y)
     */
    void dispatch_all( );
    
    //! recompose bubbles on master
    void assemble_all( );
    
    //==========================================================================
    // Fields Operations
    //==========================================================================
    
    //! MPI init ghosts exchange
    void init_exchange(); 
    
    //! MPI wait ghosts exchange
    void wait_exchange();
    
    
    //! once the bubble are dispatched, with their internal pressure
    /**
     update P and gradP.
     - P is synchronized over all the domains
     - gradP is locally computed
     */
    void compute_pressure();
    
    //! compute velocities from gradP
    /**
     - compute the local velocities
     */
    void compute_velocities();
    
    //! one function
    /**
     - dispatch_all()
     - compute_pressure
     - compute_velocities
     - synchronize U, gradP,...
     */
    void compute_fields();
    
    void advect( Real dt);
    
    Real   Lambda( Real g ) const;
    Vertex velocity_from( const Vertex &g ) const;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
    mpi::Requests   requests;    
    vector<JPack>   self_jpack;
    vector<JPack>   peer_jpack;
    void            check_borders();
};

#endif
