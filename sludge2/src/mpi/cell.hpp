#ifndef CELL_INCLUDED
#define CELL_INCLUDED

#include "../parameters.hpp"
#include "yocto/spade/mpi/workspace.hpp"
#include "yocto/spade/vtk/writer.hpp"
#include "../jpack.hpp"

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
    Array         &Bulk;      //!< count of bulk vertices per quad
    Segmenter      segmenter;
    Bubbles        bubbles;
    const Real     ymin;
    const Real     ymax;
    vtk_writer     vtk;
    vector<JPack>  jsend;
    vector<JPack>  jrecv;
    
    Real           ftol;           //!< fractional tolerance, default is 1e-5
    Real           right_pressure; //!< right pressure, default is 1
    bool           right_wall;     //!< right is a wall, default is false (set right_pressure)
    
    
    //! broadcast, find spots and bubbles, sync B field
    void dispatch( const mpi &MPI );
    
    //! save gnuplot B field
    void save_B( const string &filename ) const;
    
    //! save gnuplot B field, for outline
    void save_outB( const string &filename ) const;
    
    //! compute P gradient from bubbles + boundary conditions
    /**
     assume P is synchronized before.
     */
    void compute_gradP();
    
       
    //! compute pressure from bubbles
    /**
     - assume bubbles are pressurized by segmenter !!
     ie that the fields Penter et Pleave
     are computed.
     - compute gradP when necessary
     - compute velocities
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

    //! for debugging: VTK fields
    void save_effective(const string &filename) const;
    
    //! save where Penter.y and Pleave.y > 0, gnuplot style
    void save_effectiveY(const string &prefix) const;
    
    
    Vertex gradP_to_U( const Vertex &g ) const;
    
    //! build the bulk field
    /**
     Bulk[j][i]= #bulk vertices among B[j(+1)][i(+1)]
     */
    void build_bulk();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
    //! compute velocity
    /**
     assume grapP is computed: called at the end of compute_pressure
     */
    void compute_bulk_velocities();
    
    //! compute orhtonormal gradient
    /**
     use J->visited to do it once
     */
    void compute_junction_gn( ConstJunctionPtr J );
    void compute_spots_velocity();
    void compute_spot_velocity( Spot *spot );
};

#endif

