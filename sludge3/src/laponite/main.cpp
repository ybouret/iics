#include "visit/simulation.hpp"
#include "yocto/code/rand.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/fs/vfs.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"

int main( int argc, char *argv[] )
{
    const char *progname = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=1)
            throw exception("usage: %s config.lua", progname);
        
        Lua::State VM;
        lua_State *L = VM();
        Lua::Config::DoFile(L, argv[1]);
     
        const size_t Nx = size_t(Lua::Config::Get<lua_Number>(L, "Nx"));
        const size_t Ny = size_t(Lua::Config::Get<lua_Number>(L, "Ny"));
        
        
        
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
        
        sim.init_one_bubble("cc");
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