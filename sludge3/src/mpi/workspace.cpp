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
E1( (*this)["E1"].as<VertexArray>() ),
L1( (*this)["L1"].as<VertexArray>() ),
E2( (*this)["E2"].as<VertexArray>() ),
L2( (*this)["L2"].as<VertexArray>() ),
W( (*this)["W"].as<Array>() ),
Bulk( (*this)["Bulk"].as<Array>() ),
DeltaP( (*this)["DeltaP"].as<Array>() ),
V( (*this)["V"].as<VertexArray>() ),
right_wall(false),
P_user(0.5),
jcoll_lo(lower.y),
jcoll_up(upper.y)
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
   
    
    //==========================================================================
    // adjust collection coordinates
    //==========================================================================
    if( !MPI.IsFirst )
    {
        --(unit_t&)jcoll_lo;
    }
    
    if( !MPI.IsFinal )
    {
        ++(unit_t&)jcoll_up;
    }
}


void Workspace:: validate_bubbles(const mpi &MPI)
{
    MPI.Printf0(stderr,"\t\tvalidating....\n");
    is_valid = false;
    if( MPI.IsFirst)
    {
        bubbles.regularize_all();
        is_valid = true;
        for(const Bubble *b = bubbles.head;b;b=b->next)
        {
            // check coordinates
        }
    }
    MPI.Bcast<bool>(is_valid, 0, MPI_COMM_WORLD);
    MPI.Printf(stderr,"\t\tValidate= %s\n", is_valid ? "TRUE" : "FALSE");
}

void Workspace:: broadcast_bubbles(const mpi &MPI)
{
    MPI.Printf0(stderr, "\t\tbrodcasting...\n");
    assert(is_valid);
    ParallelBubbles::Bcast(MPI, bubbles);
}


void Workspace:: segment()
{
    junctions.load(bubbles);
    junctions.segment(B);
    
    bubbles.collect_all_markers(Y[jcoll_lo], Y[jcoll_up]);
    
    
    for( unit_t j=outline.lower.y;j<=outline.upper.y;++j)
    {
        const unit_t jm = j-1;
        const unit_t jp = j+1;
        for(unit_t i = outline.lower.x; i <= outline.upper.x; ++i)
        {
            if( B[j][i] >= 0 )
            {
                // in bubble
                Bulk[j][i]  = -1;
            }
            else
            {
                Real &bulk = Bulk[j][i];
                bulk = 0;
                if(j>outline.lower.y && B[jm][i]  < 0 ) ++bulk;
                if(j<outline.upper.y && B[jp][i]  < 0 ) ++bulk;
                if(i>outline.lower.x && B[j][i-1] < 0 ) ++bulk;
                if(i<outline.upper.x && B[j][i+1] < 0 ) ++bulk;
            }
        }
    }
    
    
}


void Workspace:: save_markers() const
{
    
    for( const Bubble *b=bubbles.head;b;b=b->next)
    {
        ios::ocstream fp( vformat("m%u.dat", unsigned(b->UID) ) , false);
        for(const Marker *m=b->markers.head;m;m=m->next)
        {
            const Tracer *tr   = m->tracer;
            const Vertex  &pos = tr->pos;
            const Vertex   tgt = pos + m->gt * tr->t; //!< tangent
            fp("%g %g\n", pos.x, pos.y);
            fp("%g %g\n", tgt.x, tgt.y);
            fp("\n");
        }
    }
    
}


bool Workspace:: initialize( const mpi &MPI )
{
    validate_bubbles(MPI);
    if(!is_valid)
    {
        return false;
    }
    
    broadcast_bubbles(MPI);
    segment();
    
    compute_pressure(MPI);
    return true;
}
