#include "visit/simulation.hpp"
#include "yocto/code/rand.hpp"
#include "yocto/string/conv.hpp"

int main( int argc, char *argv[] )
{
    try
    {
        size_t Nx = 40;
        size_t Ny = 60;
        
        if( argc > 1 )
        {
            Nx = strconv::to<size_t>( argv[1], "Nx" );
            if( argc <= 2)
                Ny = Nx;
        }
        
        if( argc > 2 )
            Ny = strconv::to<size_t>( argv[2], "Ny" );
        
        
        _rand.wseed();
        //----------------------------------------------------------------------
        //info for VisIt
        //----------------------------------------------------------------------
        const string sim_name    = "Sludge.v3";
        const string sim_comment = "Simple Laponite Using Dissolved Gas";
        const string sim_path    = ".";
        
        //----------------------------------------------------------------------
        // initialize MPI
        //----------------------------------------------------------------------
        YOCTO_MPI;
        const string       trace_name = "trace.dat";
        VisIt:: TraceFile  trace_file(MPI.CommWorldRank,trace_name );
        
        //----------------------------------------------------------------------
        // parallel VisIt
        //----------------------------------------------------------------------
        VisIt:: SetupParallel( MPI, sim_name, sim_comment, sim_path, NULL);
        
        
        //----------------------------------------------------------------------
        // Setup simulation
        //----------------------------------------------------------------------
        Simulation sim(MPI, Coord(Nx,Ny), Vertex(10,10) );
        
        __Grid::SaveDat( sim.mesh, "grid" + MPI.CommWorldID + ".dat");
        
        sim.init_one_bubble(NULL);
        sim.bubbles.gamma = 0.0;
        
        //----------------------------------------------------------------------
        // Main Loop
        //----------------------------------------------------------------------
        sim.initialize();
        VisIt:: MainLoop(sim);
        
        return 0;
    }
    catch( const exception &e )
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "Unhandled Exception" << std::endl;
    }
    return 1;
}