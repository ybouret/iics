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
U( (*this)["U"].as<VertexArray>() ),
B( (*this)["B"].as<Array>() ),
X( mesh.X() ),
Y( mesh.Y() ),
dX( mesh.dX() ),
dY( mesh.dY() ),
bubbles( sim_box ),
segmenter( mesh ),
border_segments(0),
border_peer(-1),
border_y(0),
requests( num_requests() )
{
    
    //--------------------------------------------------------------------------
    //! build the sub mesh
    //--------------------------------------------------------------------------
    mesh.regular_map_to(sim_region, sim_layout);
    bubbles.lambda = min_of( dX[dX.lower], dY[dY.lower] )/2;
    
    //--------------------------------------------------------------------------
    //! prepare the segments once the mesh if ok
    //--------------------------------------------------------------------------
    segmenter.allocate_segments();
    
    //--------------------------------------------------------------------------
    // detect parallel pbc
    //--------------------------------------------------------------------------
    if( sim_parallel )
    {
        if( 0 == sim_rank )
        {
            border_segments   = & segmenter.horizontal[lower.y];
            (int&)border_peer =  MPI.CommWorldLast;
            (Real&)border_y   = Y[lower.y];
        }
        else
        {
            if( MPI.CommWorldLast == sim_rank )
            {
                border_segments   = & segmenter.horizontal[upper.y];
                (int&)border_peer = 0; 
                (Real&)border_y   = Y[lower.y];
            }
        }
    }
    
    //! prepare ghosts
    prepare_ghosts();
    
    MPI.Printf( stderr, "Cell rank %d> region:(%7.2f->%7.2f) / layout:(%7.2f->%7.2f) / outline: (%7.2f->%7.2f) | PBC: %s (with %d)\n", sim_rank,sub_region.vmin.y, sub_region.vmax.y, Y[lower.y], Y[upper.y], Y[Y.lower], Y[Y.upper], border_segments ? "ON" : "OFF", border_peer );
    
}


void Cell:: init_exchange()
{
    _mpi::init_exchange(MPI, *this, requests);
}

void Cell:: wait_exchange()
{
    _mpi::wait_exchange(MPI, *this, requests);
}



