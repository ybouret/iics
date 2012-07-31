#include "cell.hpp"


void Cell:: check_borders()
{
    if( sim_parallel )
    {
        if( border_segments )
        {
            const int segtag = 0x5E0;
            
            //------------------------------------------------------------------
            // send/recv how many
            //------------------------------------------------------------------
            //fprintf( stderr, "%d should send %lu segments to %d\n", sim_rank, border_segments->size, border_peer );
            const size_t self_ns = border_segments->size;
            size_t       peer_ns = 0;
            MPI_Status   status;
            MPI.Sendrecv(&self_ns, sizeof(self_ns), MPI_BYTE, border_peer, segtag,
                         &peer_ns, sizeof(peer_ns), MPI_BYTE, border_peer, segtag, 
                         MPI_COMM_WORLD,status );
            //fprintf( stderr, "%d will recv %lu segments from %d\n",sim_rank,peer_ns,border_peer);
            
            //------------------------------------------------------------------
            // prepare resources
            //------------------------------------------------------------------
            JPack::HEncode(self_jpack, *border_segments); assert( self_jpack.size() == border_segments->size );
            JPack::Prepare(peer_jpack, peer_ns );         assert( peer_jpack.size() == peer_ns);
            
            //------------------------------------------------------------------
            // send resources
            //------------------------------------------------------------------
            MPI.Sendrecv(self_jpack(), self_jpack.size() * sizeof(JPack), MPI_BYTE, border_peer, segtag,
                         peer_jpack(), peer_jpack.size() * sizeof(JPack), MPI_BYTE, border_peer, segtag,
                         MPI_COMM_WORLD, status);
            
            //------------------------------------------------------------------
            // Unpack the peer border
            //------------------------------------------------------------------
            Segment::List peer_border( segmenter.s_cache );
            JPack::HDecode(peer_border, peer_jpack, segmenter, border_y, bubbles);
            
            //------------------------------------------------------------------
            // and check/merge...
            //------------------------------------------------------------------
            //MPI.__WaitFor(0.5*sim_rank);
            //fprintf( stderr, "\ncheck_borders on rank %d\n", sim_rank );
            segmenter.merge_pbc(*border_segments, border_y, peer_border, border_y);
            peer_border.empty();
        }
    }
    else 
    {
        segmenter.horizontal_pbc(lower.y, upper.y+1);
    }

}

void Cell:: dispatch_all( )
{
    MPI.Printf0( stderr, "\t---> check_and_dispatch bubbles\n");
    bubbles.check_and_dispatch_all(MPI,rescaler);
    
    MPI.Printf0( stderr, "\t---> compute bubbles properties\n");
    bubbles.check_geometries_within(Y[Y.lower], Y[Y.upper]);
    
    MPI.Printf0( stderr, "\t---> segmentation: process\n");
    segmenter.process( bubbles );
    
    MPI.Printf0( stderr, "\t---> segmentation: check_borders\n");
    check_borders();
       
    MPI.Printf0( stderr, "\t---> segmentation: assign\n");
    segmenter.assign_markers();
    
    MPI.Printf0( stderr, "\t---> segmentation: fill B\n");
    bubbles.fill(B);
}
