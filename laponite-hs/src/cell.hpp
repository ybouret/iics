#ifndef CELL_INCLUDED
#define CELL_INCLUDED 1

#include "bubbles.hpp"
#include "yocto/swamp/array2d.hpp"
#include "yocto/swamp/workspace.hpp"
#include "yocto/swamp/rmesh.hpp"

using namespace swamp;

typedef array1D<Real>         Array1D;
typedef array2D<Real>         Array;
typedef array2D<V2D>          ArrayVec;
typedef coord2D               Coord;
typedef layout<Coord>         Layout;
typedef region2D<Real>::type  Region;

typedef workspace<Layout, Real, rmesh>   WorkspaceBase;
typedef WorkspaceBase::mesh_type         Mesh;
typedef ghosts_setup<Coord>              GhostsSetup;
typedef fields_setup<Layout>             FieldsSetup;


//! base class for simulation
class Parameters : public FieldsSetup
{
public:
    //! declare all fields
    explicit Parameters(unit_t      Nx, 
                        unit_t      Ny,
                        Real        Lx,
                        Real        Ly,
                        const mpi  &MPI
                        );
    
    //! cleanup
    virtual ~Parameters() throw();
    
    const Coord  Lower;      //!< 0,0
    const Coord  Upper;      //!< Nx, Ny
    const V2D    Length;     //!< Lx,Ly
    const Layout FullLayout; //!< Lower,Upper
    const V2D    BotLeft;    //!< (0,-Ly/2)
    const V2D    TopRight;   //!< (0,Ly/2)
    const Region FullRegion; //!< BotLeft->TopRight
    const Layout SubLayout;  //!< splitted
    const Region SubRegion;  //!< spliited from SubLayout
    GhostsSetup  gs;         //!< info about ghosts
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};


#include "./intersection.hpp"



//! Hele-Shaw Cell
class Cell : public Parameters, public WorkspaceBase
{
public:
    explicit Cell(unit_t     Nx, 
                  unit_t     Ny,
                  Real       Lx,
                  Real       Ly,
                  const mpi &MPI);
    virtual ~Cell() throw();
    
    Array          &P;
    ArrayVec       &U;
    Array          &B; //!< bubble or not bubble
    const Array1D  &X;
    const Array1D  &Y;
    const Array1D  &dX;
    const Array1D  &dY;
    
    Bubbles            bubbles;
    Intersection::Pool ipool;     // cache for intersections
    Intersection::List inter;     // and intersection
    Segment::List     *horz_seg;  // along x, lower.y->upper.y
    Segment::List     *vert_seg;  // along y, lower.x->upper.x
    Segment::Pool      seg_pool;
    
    
    //! for rank 0: update_contour for each bubble
    void    master_update_topologies() throw();
    
    //! dispatch bubbles
    /**
     find their spots and their values (area,...),
     and locate each spot on the mesh.
     */
    void    dispatch_bubbles( const mpi &MPI );
    
    
    //! collect all the bubbles motions
    void    assemble_bubbles( const mpi &MPI );
    
    //! locate p->pos on the mesh
    /**
     find intersection
     */
    void locate_point( Point &p );
    
    
    //! collect points inside bubbles
    void collect_inside( vector<V2D> &pts ) const;
    
    //! advect a previously located point with field U
    void advect_point( Point &p, double dt ) const;
    
    
    //! advect all located points
    void advect_points( double dt );
    
    void save_inter( const string &filename ) const;
    void save_inside( const string &filename ) const;
    void save_grid( const string &filename ) const;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
    vector<Segment::List> segments;
    
    //! locate all points in all spots,
    void locate_points();
    void find_intersections(const V2D &P, V2D &Q, const V2D &vmin, const V2D &vmax, const U2D &pos, Bubble *bubble);
    
};



#endif