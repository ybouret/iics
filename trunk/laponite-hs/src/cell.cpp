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
    Y_SWAMP_DECL_SELF_VAR("B", Array);
    
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
B( (*this)["B"].as<Array>()    ),
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
}

void Cell:: assemble_bubbles( const mpi &MPI )
{
    bubbles.assemble_all(MPI);
}

void Cell:: save_inter( const string &filename ) const
{
    ios::ocstream fp(filename,false);
    for( const Intersection *I = inter.head; I; I=I->next )
    {
        fp("%.15g %.15g\n", I->vertex.x, I->vertex.y );
    }
}


void Cell:: save_inside( const string &filename ) const
{
    vector<V2D> pts;
    collect_inside(pts);
    ios::ocstream fp(filename,false);
    for( size_t i=pts.size();i>0;--i)
    {
        fp("%.15g %.15g\n", pts[i].x, pts[i].y);
    }
}

void Cell:: save_grid( const string &filename ) const
{
    ios::ocstream fp( filename, false );
    for( unit_t j= lower.y; j<= upper.y;++j)
    {
        if( (j&1) )
        {
            for( unit_t i= lower.x; i <= upper.x; ++i )
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
        else
        {
            for( unit_t i= upper.x; i >= lower.x; --i )
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
    }
    
    for( unit_t i=lower.x; i <= upper.x; ++i )
    {
        if( (i&1) )
        {
            for( unit_t j=lower.y; j<=upper.y;++j)
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
        else 
        {
            for( unit_t j=upper.y; j>=lower.y;--j)
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
    }
    
}



