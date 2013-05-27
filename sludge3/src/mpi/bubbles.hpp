#ifndef SLUDGE_MPI_BUBBLES_INCLUDED
#define SLUDGE_MPI_BUBBLES_INCLUDED 1

#include "../bubbles.hpp"
#include "yocto/mpi/mpi.hpp"

struct Parallel
{
    
    static void    TracerSend(const mpi &MPI, const Tracer *tr);//!< send one tracer to all slaves
    static Tracer *TracerRecv(const mpi &MPI); //!< create a tracer from master's info
    
    
    static void BubbleSend(const mpi &MPI, const Bubble *bubble); //!< send one bubble to all slaves
    static void BubbleEmit(const mpi &MPI, const Bubble *bubble); //!< emit one bubble from master
    static void BubbleRecv( const mpi &MPI, Bubbles &owner);      //!< receive one bubble from master
    
    static void BubblesEmit(const mpi &MPI,const Bubbles &bubbles); //!< send all bubbles to all slaves
    static void BubblesRecv(const mpi &MPI, Bubbles &bubbles);      //!< recv all bubbles from master
    
    
    static void BubblesBcast(const mpi &MPI, Bubbles &bubbles);
    
};


#endif

