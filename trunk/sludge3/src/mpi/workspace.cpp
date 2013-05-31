#include "workspace.hpp"
#include "yocto/code/utils.hpp"

Workspace:: ~Workspace() throw()
{
}

Workspace:: Workspace(const mpi    &MPI,
                      const Coord   N,
                      const Vertex  Q
                      ) :
Parameters(MPI,N,Q),
WorkspaceType(sim_layout,F,G),
junctions(mesh),
X( mesh.X() ),
Y( mesh.Y() ),
P( (*this)["P"].as<Array>() ),
B( (*this)["B"].as<Array>() )
{
    
    //==========================================================================
    // computing mesh
    //==========================================================================
    for(unit_t i=X.lower;i <= X.upper; ++i)
    {
        mesh.X()[i] = (delta.x)*i;
    }
    
    for(unit_t j=Y.lower;j<=Y.upper;++j)
    {
        mesh.Y()[j] = (delta.y)*j;
    }
    
    P.ld( MPI.CommWorldRank );
    
    
    //==========================================================================
    // final update
    //==========================================================================
    bubbles.lambda = min_of(delta.x,delta.y) / 2;
    
    
}


void Workspace:: broadcast_bubbles( const mpi &MPI )
{
    ParallelBubbles::Bcast(MPI, bubbles);
}

void Workspace:: segment()
{
    junctions.load(bubbles);
    junctions.segment(B);
}
