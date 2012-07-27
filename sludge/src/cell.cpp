#include "cell.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/swamp/mpi.hpp"

Cell:: ~Cell() throw()
{
}

Cell:: Cell(unit_t         Nx,
            unit_t         Ny,
            const Vertex & box,
            const mpi    & mpi_ref ) :
Parameters( Nx, Ny, box, mpi_ref),
WorkspaceBase( sub_layout, gs, *this),
P( (*this)["P"].as<Array>()    ),
B( (*this)["B"].as<Array>() ),
U( (*this)["U"].as<VertexArray>() ),
gradP( (*this)["gradP"].as<VertexArray>() ),
X( mesh.X() ),
Y( mesh.Y() ),
dX( mesh.dX() ),
dY( mesh.dY() ),
bubbles( sim_box ),
segmenter( mesh ),
delta_X( sim_box.x / Nx),
delta_Y( sim_box.y / Ny),
inv_dX2( 1/(delta_X*delta_X) ),
inv_dY2 (1/(delta_Y*delta_Y) ),
inv_two_dX( 1/(delta_X+delta_X) ),
inv_two_dY( 1/(delta_Y+delta_Y) ),
stencil_w( 1/( -2*inv_dX2 - 2*inv_dY2)),
border_segments(0),
border_peer(-1),
border_j(0),
border_y(0),
requests( num_requests() ),
self_jpack(),
peer_jpack()
{
    
    //--------------------------------------------------------------------------
    //! build the sub mesh
    //--------------------------------------------------------------------------
    {
        Array1D & _X = mesh.X();
        for(unit_t i= _X.lower; i <= _X.upper; ++i) _X[i] = i * delta_X;
        _X[upper.x] = sim_box.x;
    }
    
    {
        Array1D & _Y = mesh.Y();
        for( unit_t j=_Y.lower;j<=_Y.upper;++j) _Y[j] = bubbles.pbc.lo + j * delta_Y;
    }
    
    
    bubbles.lambda = min_of( delta_X, delta_Y )/2;
    
    //--------------------------------------------------------------------------
    //! prepare the segments once the mesh if ok
    //--------------------------------------------------------------------------
    segmenter.allocate_segments();
    
    //--------------------------------------------------------------------------
    // detect parallel pbc
    //--------------------------------------------------------------------------
    if( sim_parallel )
    {
        if( sim_master )
        {
            assert(lower.y==0);
            (unit_t&)border_j = lower.y;
            border_segments   = &segmenter.horizontal[border_j];
            (int&)border_peer = MPI.CommWorldLast;
            (Real&)border_y   = Y[border_j]; //!< bubbles.pbc.lo
        }
        else
        {
            if( sim_last )
            {
                assert(Ny==upper.y+1);
                (unit_t&)border_j = upper.y+1;
                border_segments   = & segmenter.horizontal[border_j];
                (int&)border_peer = 0;
                (Real&)border_y   = Y[border_j];
            }
        }
    }
    
    //--------------------------------------------------------------------------
    // compute the delta X/Y
    //--------------------------------------------------------------------------
    mesh.compute_deltas();
    
    
    //! prepare ghosts
    prepare_ghosts();
    
    MPI.Printf( stderr, "Cell rank %d> region:(%7.2f->%7.2f) / layout:(%7.2f->%7.2f) / outline: (%7.2f->%7.2f) | PBC: %s (with %d)\n", sim_rank,sub_region.vmin.y, sub_region.vmax.y, Y[lower.y], Y[upper.y], Y[Y.lower], Y[Y.upper], border_segments ? "ON" : "OFF", border_peer );
    
}


void Cell:: init_exchange()
{
    _mpi::init_exchange_all(MPI, *this, requests);
}

void Cell:: wait_exchange()
{
    _mpi::wait_exchange_all(MPI, *this, requests);
}




