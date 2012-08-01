#include "simulation.hpp"


int main( int argc, char *argv[] )
{
    
    try
    {
     
        AleaInit();
        
        //----------------------------------------------------------------------
        //info for VisIt
        //----------------------------------------------------------------------
        const string sim_name    = "Laponite";
        const string sim_comment = "Simple Laponite Using Dissolved Gas Experiment";
        const string sim_path    = ".";

        
        
        //----------------------------------------------------------------------
        // initialize MPI
        //----------------------------------------------------------------------
        const mpi &MPI = mpi::init( &argc, &argv );
        const int sim_rank = MPI.CommWorldRank;
        
        const string       trace_name = "trace.dat";
        VisIt:: TraceFile  trace_file(sim_rank,trace_name );
        
        //----------------------------------------------------------------------
        // parallel VisIt
        //----------------------------------------------------------------------
        VisIt:: SetupParallel( MPI, sim_name, sim_comment, sim_path, NULL);
        
        //----------------------------------------------------------------------
        // settup the simulation
        //----------------------------------------------------------------------
        Vertex box(10,20);
        //Vertex center(box.x/2,0);
        
        Simulation sim(40,70,box,MPI);
        //sim.in_walls = true;
        
        //! initialized fields and bubbles
        sim.initialize();
        
              
        //----------------------------------------------------------------------
        // Run the simulation
        //----------------------------------------------------------------------
        VisIt::MainLoop(sim);
        return 0;
    }
    catch(const exception &e)
    {
        std::cerr << e.when() << std::endl;
        std::cerr << e.what() << std::endl;
    }
    catch(...)
    {
        std::cerr << "Unhandled Exception!" << std::endl;
    }
    return  1;
}