#include "bubbles.hpp"
#include "yocto/auto-ptr.hpp"
void Parallel:: TracerSend( const mpi &MPI, const Tracer *tr )
{
    assert(tr);
    assert(0==MPI.CommWorldRank);
    for(int dest=1; dest < MPI.CommWorldSize; ++dest)
    {
        MPI.Send(&tr->pos, Tracer::NumReals, REAL_TYPE, dest, Tracer::Tag, MPI_COMM_WORLD);
    }
}



Tracer *Parallel:: TracerRecv(const yocto::mpi &MPI)
{
    assert(MPI.CommWorldRank>0);
    auto_ptr<Tracer> tr( new Tracer() );
    MPI_Status status;
    MPI.Recv(& tr->pos, Tracer::NumReals, REAL_TYPE, 0, Tracer::Tag, MPI_COMM_WORLD, status);
    return tr.yield();
}

void Parallel:: BubbleSend(const yocto::mpi &MPI, const Bubble *bubble)
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
        TracerSend(MPI, tr);
    }
}


void Parallel:: BubbleRecv( const mpi &MPI, Bubbles &owner)
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
            bubble->push_back( TracerRecv(MPI) );
        }
    }
    catch(...)
    {
        delete owner.pop_back();
        throw;
    }
    
}

static const int tag = 7;

void Parallel:: BubblesEmit(const mpi &MPI, const Bubbles &bubbles)
{
    MPI.Printf(stderr,"BubblesEmit\n");
    
    //==========================================================================
    // emit #bubbles (and extra info if needed)
    //==========================================================================
    assert(0==MPI.CommWorldRank);
    for(int dest=1;dest<MPI.CommWorldSize;++dest)
    {
        MPI.Send<size_t>(bubbles.size, dest, tag, MPI_COMM_WORLD);
    }
    return;
    
    //==========================================================================
    // emit the bubbles
    //==========================================================================
    for(const Bubble *b = bubbles.head;b;b=b->next)
    {
        BubbleSend(MPI,b);
    }
}

void Parallel:: BubblesRecv(const mpi &MPI, Bubbles &bubbles)
{
    assert(MPI.CommWorldRank>0);
    bubbles.auto_delete();
    MPI_Status status;
    MPI.Printf(stderr,"BubblesRecv\n");
    //==========================================================================
    // recv #bubbles (and extra info if needed)
    //==========================================================================
    const size_t num_bubbles = MPI.Recv<size_t>(0, tag, MPI_COMM_WORLD, status);
    MPI.Printf(stderr, "num_bubbles=%u\n",unsigned(num_bubbles) );
    return;
    
    //==========================================================================
    // recv all bubbles
    //==========================================================================
    for(size_t i=1;i<=num_bubbles;++i)
    {
        BubbleRecv(MPI, bubbles);
    }
}

void Parallel:: BubblesBcast(const mpi &MPI, Bubbles &bubbles)
{
    MPI.Printf(stderr, "BubblesBcast\n");
    if( MPI.IsFirst )
    {
        BubblesEmit(MPI,bubbles);
    }
    else
    {
        BubbleRecv(MPI,bubbles);
    }
}


