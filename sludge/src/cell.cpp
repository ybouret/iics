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
pbc_segments(0),
pbc_peer(-1),
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
            pbc_segments   = & segmenter.horizontal[lower.y];
            (int&)pbc_peer =  MPI.CommWorldLast;
        }
        else
        {
            if( MPI.CommWorldLast == sim_rank )
            {
                pbc_segments   = & segmenter.horizontal[upper.y];
                (int&)pbc_peer = 0; 
            }
        }
    }
    
    //! prepare ghosts
    prepare_ghosts();
    
    MPI.Printf( stderr, "Cell rank %d> region:(%7.2f->%7.2f) / layout:(%7.2f->%7.2f) / outline: (%7.2f->%7.2f) | PBC: %s (with %d)\n", sim_rank,sub_region.vmin.y, sub_region.vmax.y, Y[lower.y], Y[upper.y], Y[Y.lower], Y[Y.upper], pbc_segments ? "ON" : "OFF", pbc_peer );
    
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
    bubbles.check_geometries_within(Y[Y.lower], Y[Y.upper]);
    
    MPI.Printf0( stderr, "\t---> segmentation: process\n");
    segmenter.process_bubbles( bubbles );
    
    if( sim_parallel )
    {
        if( pbc_segments )
        {
            const int segtag = 0x5E0;
            fprintf( stderr, "%d should send %lu segments to %d\n", sim_rank, pbc_segments->size, pbc_peer );
            const size_t self_ns = pbc_segments->size;
            size_t       peer_ns = 0;
            MPI_Status   status;
            MPI.Sendrecv(&self_ns, sizeof(self_ns), MPI_BYTE, pbc_peer, segtag,
                         &peer_ns, sizeof(peer_ns), MPI_BYTE, pbc_peer, segtag, 
                         MPI_COMM_WORLD,status );
            fprintf( stderr, "%d will recv %lu segments from %d\n",sim_rank,peer_ns,pbc_peer);
            
            const size_t ns = self_ns + peer_ns;
            if( ns > 0 )
            {
                
            }
            
        }
    }
    else 
    {
        //segmenter.horizontal_pbc(lower.y, upper.y);
    }
    
    MPI.Printf0( stderr, "\t---> segmentation: assign\n");
    segmenter.assign_markers();
    
    MPI.Printf0( stderr, "\t---> segmentation: fill B\n");
    bubbles.fill(B);
}

void Cell:: assemble_all()
{
    bubbles.assemble_all(MPI);
}

