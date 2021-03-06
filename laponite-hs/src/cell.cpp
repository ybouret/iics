#include "cell.hpp"

Parameters:: Parameters(unit_t      Nx, 
                        unit_t      Ny,
                        Real        Lx,
                        Real        Ly,
                        const mpi & MPI) :
FieldsSetup(),
Lower(0,0),
Upper(Nx,Ny),
Length(Lx,Ly),
FullLayout( Lower, Upper),
BotLeft(0,-Length.y/2),
TopRight(Length.x,Length.y/2),
FullRegion(BotLeft,TopRight),
SubLayout( FullLayout.split(MPI.CommWorldRank, MPI.CommWorldSize) ),
SubRegion( FullRegion.split(MPI.CommWorldRank, MPI.CommWorldSize) ),
gs()
{
    assert(Lx>0);
    assert(Ly>0);
    Y_SWAMP_DECL_SELF_VAR("P",    Array);
    Y_SWAMP_DECL_SELF_VAR("U",    ArrayVec);
    Y_SWAMP_DECL_SELF_AUX("B",    Array);

    if( MPI.IsParallel )
    {
        gs.lower.count.y = 2;
        gs.lower.peer.y  = MPI.CommWorldPrev();
        gs.upper.count.y = 2;
        gs.upper.peer.y  = MPI.CommWorldNext();
    }
    else 
    {
        gs.local.count.y = 2;
    }
}

Parameters:: ~Parameters() throw()
{
}


Cell:: Cell(unit_t      Nx, 
            unit_t      Ny,
            Real        Lx,
            Real        Ly,
            const mpi & MPI ) :
Parameters( Nx, Ny, Lx, Ly, MPI),
WorkspaceBase( SubLayout, gs, *this),
P( (*this)["P"].as<Array>()    ),
U( (*this)["U"].as<ArrayVec>() ),
B( (*this)["B"].as<Array>() ),
X( mesh.X() ),
Y( mesh.Y() ),
dX( mesh.dX() ),
dY( mesh.dY() ),
bubbles( Length ),
ipool(),
inter(ipool),
horz_seg(0),
vert_seg(0),
seg_pool(),
segments()
{
    //! build the sub mesh
    mesh.regular_map_to(FullRegion, FullLayout);
    const V2D v( dX[lower.x], dY[lower.y]);
    (Real &)(bubbles.lambda) = 0.5 * v.norm();
    
    //! prepare the segments: for each X/Y
    const Segment::List seg( seg_pool );
    segments.make(X.width+Y.width ,seg);
    horz_seg = &segments[1] - Y.lower;
    vert_seg = &segments[1+Y.width] - X.lower;
}

Cell:: ~Cell() throw()
{
    
}

void Cell:: master_update_topologies() throw()
{
    bubbles.update_topologies();
}


void Cell:: dispatch_bubbles( const mpi &MPI ) 
{
    //! broadcast
    MPI.Printf0(stderr, "\t*** dispatch all bubbles\n");
    bubbles.dispatch_all(MPI);
    
    //! compute spots and local properties
    MPI.Printf0(stderr, "\t*** build spots and values\n");
    bubbles.spots_and_values_within(SubRegion.vmin.y, SubRegion.vmax.y);
    
    //! locate points
    MPI.Printf0(stderr, "\t*** locate points\n");
    locate_points();
    
    //! locate markers
    MPI.Printf0(stderr, "\t*** locate inside markers\n");
    locate_markers();
    
    //! propagate and build B
    MPI.Printf0(stderr, "\t*** propagate inside markers\n");
    propagate_markers(MPI);
    
}

void Cell:: assemble_bubbles( const mpi &MPI )
{
    bubbles.assemble_all(MPI);
}




