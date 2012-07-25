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
requests( num_requests() )
{
    
    //! build the sub mesh
    mesh.regular_map_to(sim_region, sim_layout);
    bubbles.lambda = min_of( dX[dX.lower], dY[dY.lower] )/2;
    
    //! prepare the segments once the mesh if ok
    segmenter.allocate_segments();
    
    //! prepare ghosts
    prepare_ghosts();
    
    MPI.Printf( stderr, "Cell rank %d> \n", sim_rank);
    
}


void Cell:: init_exchange()
{
    _mpi::init_exchange(MPI, *this, requests);
}

void Cell:: wait_exchange()
{
    _mpi::wait_exchange(MPI, *this, requests);
}


void Cell:: dispatch_all( )
{
    MPI.Printf0( stderr, "\t---> check_and_dispatch bubbles\n");
    bubbles.check_and_dispatch_all(MPI);
    
    MPI.Printf0( stderr, "\t---> compute bubbles properties\n");
    bubbles.check_geometries_within(sub_region.vmin.y, sub_region.vmax.y);
    
    
}

void Cell:: assemble_all()
{
    bubbles.assemble_all(MPI);
}

