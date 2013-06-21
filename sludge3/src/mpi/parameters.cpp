#include "parameters.hpp"

Parameters:: ~Parameters() throw() {}

Parameters:: Parameters(const mpi    &MPI,
                        const Coord  &N,
                        const Vertex &Q) :
F( 32 ),
G(),
full_layout( Coord(0,0), N),
sim_layout( full_layout.split(MPI.CommWorldRank, MPI.CommWorldSize,1) ),

delta(     Q.x/N.x,             Q.y/N.y        ),
two_delta( delta.x+delta.x,     delta.y+delta.y),
inv_delta( 1/delta.x,           1/delta.y      ),
order1fac( 1/two_delta.x,       1/two_delta.y  ),
order2fac( 1/(delta.x*delta.x), 1/(delta.y*delta.y) ),

weight( -2*order2fac.y, -2*order2fac.y),
rb_weight( weight.x + weight.y ),

full_region( Vertex(0,0), Q ),
sim_region( Vertex(0, sim_layout.lower.y * delta.y), Vertex(Q.x,sim_layout.upper.y*delta.y) ),
bulk_imin( sim_layout.lower.x+1),
bulk_imax( sim_layout.upper.x-1),
bulk_jmin( sim_layout.lower.y + (MPI.IsFirst ? 1 : 0 ) ),
bulk_jmax( sim_layout.upper.y - (MPI.IsFinal ? 1 : 0 ) ),

ftol(1e-5)
{
    
    MPI.Printf(stderr, "Full Layout: (%ld,%ld) -> (%ld,%ld)\n", full_layout.lower.x, full_layout.lower.y, full_layout.upper.x, full_layout.upper.y);
    MPI.Printf(stderr, "Sim  Layout: (%ld,%ld) -> (%ld,%ld)\n", sim_layout.lower.x, sim_layout.lower.y, sim_layout.upper.x, sim_layout.upper.y);
    MPI.Printf(stderr, "Sim  Region: [%g %g]' --> [%g %g]'\n", sim_region.vmin.x, sim_region.vmin.y, sim_region.vmax.x, sim_region.vmax.y);
    MPI.Printf(stderr, "Red/Black Factor= %g\n", 1/rb_weight);
    
    //==========================================================================
    // prepare fields
    //==========================================================================
    Y_SPADE_FIELD(F, "P",      Array);
    Y_SPADE_FIELD(F, "B",      Array);
    Y_SPADE_FIELD(F, "gradP",  VertexArray);
    Y_SPADE_FIELD(F, "E1",     VertexArray);
    Y_SPADE_FIELD(F, "L1",     VertexArray);
    Y_SPADE_FIELD(F, "E2",     VertexArray);
    Y_SPADE_FIELD(F, "L2",     VertexArray);
    Y_SPADE_FIELD(F, "DeltaP", Array);
    Y_SPADE_FIELD(F, "W",      Array);
    Y_SPADE_FIELD(F, "Bulk",   Array);
    Y_SPADE_FIELD(F, "V",      VertexArray);
    
    //==========================================================================
    // MPI connectivity
    //==========================================================================
    if( MPI.IsParallel )
    {
        if( MPI.IsFirst )
        {
            G.set_async(ghost::at_upper_y, NumGhosts, MPI.CommWorldNext());
        }
        else
        {
            if(MPI.IsFinal)
            {
                G.set_async(ghost::at_lower_y, NumGhosts, MPI.CommWorldPrev());
            }
            else
            {
                G.set_async(ghost::at_upper_y, NumGhosts, MPI.CommWorldNext());
                G.set_async(ghost::at_lower_y, NumGhosts, MPI.CommWorldPrev());
            }
        }
    }
    else
    {
        // no ghosts !
    }
    
    //==========================================================================
    // Markers Collection
    //==========================================================================
    
}