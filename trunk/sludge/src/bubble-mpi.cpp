#include "bubble.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
void Bubble:: dispatch_topology( const mpi &MPI )
{
    if( MPI.CommWorldSize >0 )
    {
        const bool master     = MPI.IsMaster;
        size_t     num_tracer = 0;
        
        //--------------------------------------------------------------------------
        // initialize
        //--------------------------------------------------------------------------
        if( master )
        {
            num_tracer = size;
        }
        else 
        {
            empty();
        }
        
        //--------------------------------------------------------------------------
        // broadcast #tracers
        //--------------------------------------------------------------------------
        MPI.BcastAs<size_t>(num_tracer, 0, MPI_COMM_WORLD);
        MPI.Bcast( &pressure, 1, MPI_REAL_TYPE, 0, MPI_COMM_WORLD);
        MPI.Printf(stderr, "rank %d> #num_tracer=%lu\n", MPI.CommWorldRank, num_tracer);
        if(!master)
            append(num_tracer);
        
        //--------------------------------------------------------------------------
        // broadcast topology
        //--------------------------------------------------------------------------
        Tracer *p = root;
        for( size_t i=num_tracer;i>0;--i,p=p->next)
        {
            MPI.Bcast( &p->vertex, Tracer::IO_COUNT, MPI_REAL_TYPE, 0, MPI_COMM_WORLD);
        }
    }
    
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
void Bubble:: assemble_topology( const mpi &MPI )
{
    static const int tag = 1000;
    if( MPI.CommWorldSize >0 )
    {
        const bool master   = MPI.IsMaster;
        //MPI.Barrier(MPI_COMM_WORLD);
        if( master )
        {
            MPI_Status status;
            fprintf( stderr, "master changed=%u\n", unsigned(spots.size) );
            for( int r = 1; r < MPI.CommWorldSize; ++r )
            {
                //--------------------------------------------------------------
                // receive the num_changed
                //--------------------------------------------------------------
                
                const size_t num_changed = MPI.RecvAs<size_t>(r, tag, MPI_COMM_WORLD, status);
                fprintf( stderr, "rank %d changed=%u\n", r, unsigned(num_changed) );
                
                //--------------------------------------------------------------
                // far all changed
                //--------------------------------------------------------------
                Tracer *p = root;
                for( size_t i=0; i < num_changed; ++i )
                {
                    //----------------------------------------------------------
                    // receive the #jump to perform
                    //----------------------------------------------------------
                    size_t jump = MPI.RecvAs<size_t>(r, tag, MPI_COMM_WORLD, status);
                    
                    //----------------------------------------------------------
                    // move forward
                    //----------------------------------------------------------
                    while(jump-->0) 
                    { 
                        assert(p!=NULL);
                        p=p->next; 
                    }
                    assert(p!=NULL);
                    
                    //----------------------------------------------------------
                    // receive the new coordinates
                    //----------------------------------------------------------
                    MPI.Recv(& p->vertex, 2, MPI_REAL_TYPE, r, tag, MPI_COMM_WORLD, status);
                }
            }
        }
        else 
        {
            //------------------------------------------------------------------
            // send the num_changed
            //------------------------------------------------------------------
            const size_t num_changed = spots.size;
            MPI.SendAs<size_t>(num_changed, 0, tag, MPI_COMM_WORLD);
            
            //------------------------------------------------------------------
            // far all changed
            //------------------------------------------------------------------
            const Spot *spot = spots.head;
            for( size_t i=0; i<num_changed; ++i , spot=spot->next)
            {
                //--------------------------------------------------------------
                // send the #jump to perform
                //--------------------------------------------------------------
                MPI.SendAs<size_t>(spot->jump, 0, tag, MPI_COMM_WORLD);
                
                //--------------------------------------------------------------
                // send the new coordinates
                //--------------------------------------------------------------
                MPI.Send(& spot->handle->vertex, 2, MPI_REAL_TYPE, 0, tag, MPI_COMM_WORLD);
            }
        }
    }
    
}

#if 0
////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
static const int __marker_tag = 666;

static inline void __send_marker( const Marker *g, int dest, const mpi &MPI )
{
    MPI.Send( & g->coord, sizeof(Coord), MPI_BYTE, dest, __marker_tag, MPI_COMM_WORLD);
}

static inline void __recv_marker(  Marker *g, int from, const mpi &MPI )
{
    MPI_Status status;
    MPI.Recv( & g->coord, sizeof(Coord), MPI_BYTE, from, __marker_tag, MPI_COMM_WORLD, status);
}


static inline void __exch_marker(Bubble    &bubble, 
                                 const int  i_send,
                                 const int  i_recv, 
                                 const mpi &MPI )
{
    MPI_Status  status;
    
    //--------------------------------------------------------------------------
    // send how many I will send
    //--------------------------------------------------------------------------
    const Marker::List  &self   = bubble.markers;
    const size_t         n_send = self.size;
    MPI.Printf(stderr, "rank %d> send #=%lu to %d\n", MPI.CommWorldRank, n_send, i_send);
    MPI.SendAs<size_t>(n_send, 
                       i_send, 
                       __marker_tag, 
                       MPI_COMM_WORLD);
    
    //--------------------------------------------------------------------------
    // recv how many I will receive
    //--------------------------------------------------------------------------
    const size_t n_recv = MPI.RecvAs<size_t>(i_recv, 
                                             __marker_tag, 
                                             MPI_COMM_WORLD, 
                                             status);
    MPI.Printf( stderr, "rank %d> recv #=%lu from %d\n", MPI.CommWorldRank, n_recv, i_recv);
    
    //--------------------------------------------------------------------------
    // send own markers
    //--------------------------------------------------------------------------
    for( const Marker *g = self.head; g; g=g->next )
    {
        __send_marker(g, i_send, MPI);
    }
    
    //--------------------------------------------------------------------------
    // collect neighbor's markers
    //--------------------------------------------------------------------------
    Marker::List &peer = bubble.borders;
    for( size_t i=n_recv; i>0; --i )
    {
        __recv_marker( peer.append(), i_recv, MPI);
    }
}

void Bubble:: propagate_markers(const yocto::mpi &MPI)
{
    if( MPI.IsParallel )
    {
        borders.empty();
        const int i_prev = MPI.CommWorldPrev();
        const int i_next = MPI.CommWorldNext();            
        __exch_marker(*this,i_prev,i_next,MPI);
        if(i_prev!=i_next)
            __exch_marker(*this, i_next, i_prev, MPI);
        markers.merge_back(borders);
    }
}
#endif
