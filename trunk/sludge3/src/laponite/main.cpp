#include "visit/simulation.hpp"

int main( int argc, char *argv[] )
{
    try
    {
        //--------------------------------------------------------------------------
        //info for VisIt
        //--------------------------------------------------------------------------
        const string sim_name    = "Sludge.v3";
        const string sim_comment = "Simple Laponite Using Dissolved Gas";
        const string sim_path    = ".";
        
        //--------------------------------------------------------------------------
        // initialize MPI
        //--------------------------------------------------------------------------
        YOCTO_MPI;
        const string       trace_name = "trace.dat";
        VisIt:: TraceFile  trace_file(MPI.CommWorldRank,trace_name );
        
        //--------------------------------------------------------------------------
        // parallel VisIt
        //--------------------------------------------------------------------------
        VisIt:: SetupParallel( MPI, sim_name, sim_comment, sim_path, NULL);

        
        Simulation sim(MPI, Coord(20,30), Vertex(10,10) );
        
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