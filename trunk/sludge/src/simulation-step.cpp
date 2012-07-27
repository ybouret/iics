#include "simulation.hpp"

void Simulation:: step()
{
    VisIt::Simulation::step();
    
    // move concerned points
    bubbles.first()->translate(Vertex(0,0.2*delta_Y));
    advect(0.1);
    
    // send back information to master
    assemble_all();
    
    
    // - process topologies on master
    // - update B, compute P, gradP, U
    // - and all exchanges
    compute_fields();
    
    
}

