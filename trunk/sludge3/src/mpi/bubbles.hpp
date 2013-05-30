#ifndef SLUDGE_MPI_BUBBLES_INCLUDED
#define SLUDGE_MPI_BUBBLES_INCLUDED 1

#include "../bubbles.hpp"
#include "yocto/mpi/mpi.hpp"

class ParallelBubbles
{
public:
    // send bubble to slaves
    static void Bcast(const mpi &MPI, Bubbles &bubbles);
    
    // gather info from slaves
    static void Collect(const mpi &MPI, Bubbles &bubles);
    
private:
    static void Send(const mpi &MPI,const Bubbles &bubbles); //!< send all bubbles to all slaves
    static void Recv(const mpi &MPI, Bubbles &bubbles);      //!< recv all bubbles from master
};

struct ParallelTracer
{
    static void    Send(const mpi &MPI, const Tracer *tr); //!< send one tracer to all slaves
    static Tracer *Recv(const mpi &MPI);                   //!< create a tracer from master's info
};

struct ParallelBubble
{
    static void Send(const mpi &MPI, const Bubble *bubble); //!< send one bubble to all slaves
    static void Recv( const mpi &MPI, Bubbles &owner);      //!< receive one bubble from master
    
    static void SendMarkers( const mpi &MPI, const Bubble *bubble); //!< send one bubble markers to master
    static void RecvMarkers( const mpi &MPI, Bubble *bubble);       //!< recv markers from all bubbles
};

#endif

