#include "bubbles.hpp"
#include "yocto/auto-ptr.hpp"

////////////////////////////////////////////////////////////////////////////////
//
// PARALLEL TRACER
//
////////////////////////////////////////////////////////////////////////////////
void ParallelTracer:: Send( const mpi &MPI, const Tracer *tr )
{
    assert(tr);
    assert(0==MPI.CommWorldRank);
    for(int dest=1; dest < MPI.CommWorldSize; ++dest)
    {
        MPI.Send(&tr->pos, Tracer::NumReals, REAL_TYPE, dest, Tracer::Tag, MPI_COMM_WORLD);
    }
}



Tracer *ParallelTracer:: Recv(const yocto::mpi &MPI)
{
    assert(MPI.CommWorldRank>0);
    auto_ptr<Tracer> tr( new Tracer() );
    MPI_Status status;
    MPI.Recv(& tr->pos, Tracer::NumReals, REAL_TYPE, 0, Tracer::Tag, MPI_COMM_WORLD, status);
    return tr.yield();
}


////////////////////////////////////////////////////////////////////////////////
//
// PARALLEL BUBBLE
//
////////////////////////////////////////////////////////////////////////////////
void ParallelBubble:: Send(const yocto::mpi &MPI, const Bubble *bubble)
{
    assert(bubble);
    assert(0==MPI.CommWorldRank);
    assert(bubble->size>=3);
    
    //==========================================================================
    //-- send how many tracers and extra info
    //==========================================================================
    for(int dest=1;dest<MPI.CommWorldSize;++dest)
    {
        MPI.Send(bubble->size, dest, Bubble::Tag, MPI_COMM_WORLD);
        MPI.Send(&bubble->G, Bubble::NumReals, REAL_TYPE, dest, Bubble::Tag, MPI_COMM_WORLD);
    }
    
    //==========================================================================
    //-- send individual tracers
    //==========================================================================
    const Tracer *tr = bubble->root;
    for(size_t i=bubble->size;i>0;--i,tr=tr->next)
    {
        ParallelTracer:: Send(MPI, tr);
    }
}


void ParallelBubble:: Recv( const mpi &MPI, Bubbles &owner)
{
    assert(MPI.CommWorldRank>0);
    
    //==========================================================================
    // create the bubble
    //==========================================================================
    Bubble *bubble = owner.append();
    try
    {
        MPI_Status status;
        //======================================================================
        // receive the #tracers
        //======================================================================
        const size_t num_tracers = MPI.Recv<size_t>(0, Bubble::Tag, MPI_COMM_WORLD, status);
        if( num_tracers < 3)
            throw exception("BubbleRecv(#tracers<3)");
        
        //======================================================================
        // receive the extra data
        //======================================================================
        MPI.Recv(& bubble->G, Bubble::NumReals, REAL_TYPE, 0, Bubble::Tag, MPI_COMM_WORLD, status);
        
        //======================================================================
        // receive the tracers
        //======================================================================
        for(size_t i=1;i<=num_tracers;++i)
        {
            bubble->push_back( ParallelTracer:: Recv(MPI) );
        }
    }
    catch(...)
    {
        delete owner.pop_back();
        throw;
    }
    
}

void ParallelBubble:: SendMarkers(const yocto::mpi &MPI, const Bubble *bubble)
{
    assert(bubble);
    assert(bubble->size>=3);
    assert(MPI.CommWorldRank>0);
    
    //==========================================================================
    // send #markers
    //==========================================================================
    assert(bubble->markers.size<=bubble->size);
    MPI.Send<size_t>(bubble->markers.size, 0, Marker::Tag, MPI_COMM_WORLD);
    
    //==========================================================================
    // send all markers
    //==========================================================================
    for( const Marker *mk=bubble->markers.head;mk;mk=mk->next)
    {
        //----------------------------------------------------------------------
        // send shift
        //----------------------------------------------------------------------
        MPI.Send<size_t>(mk->shift,0,Marker::Tag, MPI_COMM_WORLD);
        
        //----------------------------------------------------------------------
        // send info
        //----------------------------------------------------------------------
        MPI.Send(& mk->tracer->pos, 2, REAL_TYPE, 0, Marker::Tag, MPI_COMM_WORLD);
    }
    
}

void ParallelBubble:: RecvMarkers(const yocto::mpi &MPI, Bubble *bubble)
{
    //==========================================================================
    // recv #markers
    //==========================================================================
    assert(bubble);
    assert(bubble->size>0);
    assert(0==MPI.CommWorldRank);
    MPI_Status status;
    for( int source=1; source < MPI.CommWorldSize; ++source )
    {
        const size_t m = MPI.Recv<size_t>(source, Marker::Tag, MPI_COMM_WORLD, status);
        if(m>bubble->size)
            throw exception("RecvMarkers(TOO MANY)");
        Tracer *tr = bubble->root;
        for(size_t i=1; i<=m; ++i )
        {
            //------------------------------------------------------------------
            // recv shift
            //------------------------------------------------------------------
            size_t shift = MPI.Recv<size_t>(source, Marker::Tag, MPI_COMM_WORLD, status);
            while(shift-->0)
            {
                tr=tr->next;
            }
            
            //------------------------------------------------------------------
            // recv info
            //------------------------------------------------------------------
            MPI.Recv(& tr->pos, 2, REAL_TYPE, source, Marker::Tag, MPI_COMM_WORLD, status);
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
//
// PARALLEL BUBBLES
//
////////////////////////////////////////////////////////////////////////////////
void ParallelBubbles:: Send(const mpi &MPI, const Bubbles &bubbles)
{
    
    //==========================================================================
    // emit #bubbles (and extra info if needed)
    //==========================================================================
    assert(0==MPI.CommWorldRank);
    for(int dest=1;dest<MPI.CommWorldSize;++dest)
    {
        MPI.Send<size_t>(bubbles.size, dest, Bubbles::Tag, MPI_COMM_WORLD);
    }
    
    //==========================================================================
    // emit the bubbles
    //==========================================================================
    for(const Bubble *b = bubbles.head;b;b=b->next)
    {
        ParallelBubble:: Send(MPI,b);
    }
}

void ParallelBubbles:: Recv(const mpi &MPI, Bubbles &bubbles)
{
    assert(MPI.CommWorldRank>0);
    bubbles.auto_delete();
    MPI_Status status;
    //==========================================================================
    // recv #bubbles (and extra info if needed)
    //==========================================================================
    const size_t num_bubbles = MPI.Recv<size_t>(0, Bubbles::Tag, MPI_COMM_WORLD, status);
    
    //==========================================================================
    // recv all bubbles
    //==========================================================================
    for(size_t i=1;i<=num_bubbles;++i)
    {
        ParallelBubble:: Recv(MPI, bubbles);
    }
}

void ParallelBubbles:: Bcast(const mpi &MPI, Bubbles &bubbles)
{
    if( MPI.IsFirst )
    {
        ParallelBubbles:: Send(MPI,bubbles);
    }
    else
    {
        ParallelBubbles:: Recv(MPI,bubbles);
    }
}


void ParallelBubbles:: Collect(const mpi &MPI, Bubbles &bubbles)
{
    
    if( MPI.IsFirst )
    {
        for( Bubble *b = bubbles.head;b;b=b->next)
        {
            ParallelBubble::RecvMarkers(MPI,b);
        }
    }
    else
    {
        for(const Bubble *b=bubbles.head;b;b=b->next)
        {
            ParallelBubble::SendMarkers(MPI,b);
        }
    }

}


