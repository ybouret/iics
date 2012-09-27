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
        const Coord  N(20,40);
        const Vertex L(3.0,5.0);
        Simulation   sim(MPI,N,L);
        MPI.PrintfI(stderr, "layout: (%d,%d) -> (%d,%d)\n", int(sim.lower.x),int(sim.lower.y),int(sim.upper.x),int(sim.upper.y));
        MPI.__WaitFor(0.1);
        for( unsigned i=1; i <= sim.async_count; ++i )
        {
            MPI.PrintfI(stderr, "\t\tAsyncGhost #%u: peer=%d\n",i, sim.get_async(i).peer);
            MPI.__WaitFor(0.1);
        }
        
        //----------------------------------------------------------------------
        // First time init
        //----------------------------------------------------------------------
        sim.initialize();
        
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