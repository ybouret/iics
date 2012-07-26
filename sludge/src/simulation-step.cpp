#include "simulation.hpp"

void Simulation:: step()
{
    VisIt::Simulation::step();
    
    // move concerned points
    //advect_points(1);
    bubbles.first()->translate(Vertex(0,0.1*Alea()));
    
    // send back information to master
    assemble_all();
    
    
    //process topologies on master
    dispatch_all();
    compute_pressure();
    
#if 0
    init_exchange();
    wait_exchange();
#endif
    
}

