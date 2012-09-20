#include "cell.hpp"

Cell:: ~Cell() throw()
{
    
}

Cell:: Cell(const mpi         &MPI,
            const Coord       &N,
            const Vertex      &L,
            const FieldsSetup &F ) :
Parameters(N,L,MPI.CommWorldRank,MPI.CommWorldSize),
Workspace( sim_layout, F, sim_ghosts ),
segmenter( mesh ),
bubbles( pbc )
{
    //-- compute the grid
    setup_grid( mesh );
    
    //-- create segmenter
    segmenter.create();
    
}

