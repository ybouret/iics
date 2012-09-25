#include "../visit/simulation.hpp"

int main( int argc, char *argv[] )
{
    try
    {
        //----------------------------------------------------------------------
        //info for VisIt
        //----------------------------------------------------------------------
        const string sim_name    = "Laponite";
        const string sim_comment = "Laponite w/ Dissolved Gaz";
        const string sim_path    = ".";
        
        //----------------------------------------------------------------------
        // Setup MPI
        //----------------------------------------------------------------------
        YOCTO_MPI;
        
        //----------------------------------------------------------------------
        // Trace File
        //----------------------------------------------------------------------
        const string       trace_name = "trace.dat";
        VisIt:: TraceFile  trace_file(MPI.CommWorldRank,trace_name);
        
        //----------------------------------------------------------------------
        // Setup VisIt
        //----------------------------------------------------------------------
        const string *sim_interface = 0;
        VisIt:: SetupParallel( MPI, sim_name, sim_comment, sim_path, sim_interface);
        
        
        //----------------------------------------------------------------------
        // Setup Simulation
        //----------------------------------------------------------------------
        Simulation sim(MPI);
        
        //----------------------------------------------------------------------
        // Process VisIt main loop
        //----------------------------------------------------------------------
        VisIt::MainLoop(sim);
        
        return 0;
    }
    catch( const exception &e )
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "unhandled exception!" << std::endl;
    }
    return 1;
}