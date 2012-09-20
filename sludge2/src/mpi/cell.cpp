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
X( mesh.X() ),
Y( mesh.Y() ),
segmenter( mesh ),
bubbles( pbc ),
ymin(0),
ymax(0)
{
    
    //-- compute the grid
    setup_grid( mesh );
    
    //-- tune the data
    (Real&)ymin = Y[ this->lower.y ];
    (Real&)ymax = Y[ this->upper.y ];
    MPI.PrintfI(stderr, "ymin=%g / ymax=%g / pbc=[%g;%g]\n", ymin, ymax, pbc.lo, pbc.up);
    
    //-- create segments
    segmenter.create();
    
    
}

void Cell:: dispatch( const mpi &MPI )
{
    MPI.Printf0(stderr, "\tdispatch bubbles...\n");
    bubbles.dispatch(MPI);
    MPI.Printf0(stderr, "\tlocate spots...\n");
    for( Bubble *b = bubbles.first(); b;b=b->next)
    {
        b->locate_spots(ymin, ymax);
    }
}
