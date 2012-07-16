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
    Y_SWAMP_DECL_SELF_VAR("P", Array);
    Y_SWAMP_DECL_SELF_VAR("U", ArrayVec);
    
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
bubbles( Length.y )
{
    //! build the sub mesh
    mesh.regular_map_to(FullRegion, SubLayout);
        
}

Cell:: ~Cell() throw()
{
    
}

void Cell:: update_topologies() throw()
{
    bubbles.update_all_points();
}


void Cell:: dispatch_bubbles( const mpi &MPI ) 
{
    bubbles.dispatch_all(MPI);
    bubbles.update_all_spots(SubRegion.vmin.y, SubRegion.vmax.y);
}

void Cell:: collect_bubbles( const mpi &MPI )
{
    bubbles.collect_all(MPI);
}
