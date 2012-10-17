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
        const Coord  N(30,40);
        const Vertex L(6.0,5.0);
        Simulation   sim(MPI,N,L);
        MPI.PrintfI(stderr, "layout: (%d,%d) -> (%d,%d)\n", int(sim.lower.x),int(sim.lower.y),int(sim.upper.x),int(sim.upper.y));
        MPI.WaitFor(0.1);
      
        SaveGrid( sim.mesh, vformat("g%d.%d.dat",sim.par_size,sim.par_rank));
        
        sim.bubbles.gamma = 0.1;
        
        //----------------------------------------------------------------------
        // First time init
        //----------------------------------------------------------------------
        sim.initialize();
        
        //sim.segmenter.save_vtk_n( "j.vtk");
        //sim.bubbles.first()->save_vtk("b.vtk");
        //sim.bubbles.first()->save_vtk_n("b-n.vtk");
        
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