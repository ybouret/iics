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
    
    Array    &P;
    ArrayVec &U;
    Array1D  &X;
    Array1D  &Y;
    
    
    Bubbles bubbles;
    
    //! for rank 0: update_contour for each bubble
    void    master_update_topologies() throw();
    
    //! dispatch bubbles, find their spots and their values (area,...)
    void    dispatch_bubbles( const mpi &MPI );
    
    //! collect all the bubbles motions
    void    assemble_bubbles( const mpi &MPI );
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Cell);
};



#endif
