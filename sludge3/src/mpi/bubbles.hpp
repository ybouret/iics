#ifndef SLUDGE_MPI_BUBBLES_INCLUDED
#define SLUDGE_MPI_BUBBLES_INCLUDED 1

#include "../bubbles.hpp"
#include "yocto/mpi/mpi.hpp"

struct Parallel
{
    //! send one tracer to all slaves
    /**
     */
    void    TracerSend(const mpi &MPI, const Tracer *tr);
    
    
    //! create a tracer from master's info
    /**
     assuming rank > 0
     */
    Tracer *TracerRecv(const mpi &MPI);
    
    
    //! send one bubble to all slaves
    /**
     \param dest > 0
     assuming rank == 0
    */
    void BubbleSend(const mpi &MPI, const Bubble *bubble);
    
    //! send one bubble to all slaves
    /**
     assuming rank = 0;
     */
    void BubbleEmit(const mpi &MPI, const Bubble *bubble);
    
    //! receive one bubble from master
    /** 
     assuming rank>0
     */
    void BubbleRecv( const mpi &MPI, Bubbles &owner);
    
    
    //! send all bubbles to all slaves
    void BubblesEmit(const mpi &MPI,const Bubbles &bubbles);
    
    //! recv all bubbles from master
    void BubblesRecv(const mpi &MPI, Bubbles &bubbles);
    
    
};


#endif

