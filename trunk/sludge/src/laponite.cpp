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
        Vertex box(10,10);
        Vertex center(box.x/2,0);
        
        Simulation sim(20,30,box,MPI);
        
        //! initialized fields and bubbles
        sim.initialize();
        
        //! dispatch bubbles/geometry
        sim.dispatch_all();
        
        //! prepare fields
        sim.init_exchange();
        sim.wait_exchange();
        
        //----------------------------------------------------------------------
        // Run the simulation
        //----------------------------------------------------------------------
        VisIt::MainLoop(sim);
        return 0;
    }
    catch(const exception &e)
    {
        
    }
    catch(...)
    {
        std::cerr << "Unhandled Exception!" << std::endl;
    }
    return  1;
}