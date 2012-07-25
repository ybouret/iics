#include "parameters.hpp"

Parameters:: ~Parameters() throw()
{
}

Parameters:: Parameters(unit_t        Nx, 
                        unit_t        Ny,
                        const Vertex &box,
                        const mpi    &mpi_ref) :
FieldsSetup(),
MPI(mpi_ref),
sim_rank(     MPI.CommWorldRank ),
sim_size(     MPI.CommWorldSize ),
sim_master(   MPI.IsMaster      ),
sim_parallel( MPI.IsParallel    ),
sim_lower(0,0),
sim_upper(Nx,Ny),
sim_box(box),
sim_layout( sim_lower, sim_upper),
sim_lower_corner(0,-box.y/2),
sim_upper_corner(box.x,box.y/2),
sim_region(sim_lower_corner,sim_upper_corner),
sub_layout(sim_layout.split(sim_rank, sim_size) ),
sub_region(sim_region.split(sim_rank, sim_size) ),
gs()
{
    assert(sim_box.x>0);
    assert(sim_box.y>0);
    Y_SWAMP_DECL_SELF_VAR("P",    Array);
    Y_SWAMP_DECL_SELF_VAR("U",    VertexArray);
    Y_SWAMP_DECL_SELF_AUX("B",    Array);
    
    if( sim_parallel )
    {
        gs.lower.count.y = 2;
        gs.lower.peer.y  = MPI.CommWorldPrev();
        gs.upper.count.y = 2;
        gs.upper.peer.y  = MPI.CommWorldNext();
    }
    else 
    {
        gs.local.count.y = 2;
    }
}
