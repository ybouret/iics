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
is_valid(true),
X( mesh.X() ),
Y( mesh.Y() ),
P( (*this)["P"].as<Array>() ),
B( (*this)["B"].as<Array>() ),
gradP( (*this)["gradP"].as<VertexArray>() ),
Enter( (*this)["Enter"].as<VertexArray>() ),
Leave( (*this)["Leave"].as<VertexArray>() )
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


void Workspace:: validate_bubbles(const mpi &MPI)
{
    
    is_valid = false;
    if( MPI.IsFirst)
    {
        bubbles.regularize_all();
        is_valid = true;
        for(const Bubble *b = bubbles.head;b;b=b->next)
        {
            
        }
        MPI.Bcast(is_valid, 0, MPI_COMM_WORLD);
    }
    MPI.Printf(stderr,"Validate= %s\n", is_valid ? "TRUE" : "FALSE");
}

void Workspace:: broadcast_bubbles(const mpi &MPI)
{
    assert(is_valid);
    ParallelBubbles::Bcast(MPI, bubbles);
}


void Workspace:: segment()
{
    junctions.load(bubbles);
    junctions.segment(B);
}
