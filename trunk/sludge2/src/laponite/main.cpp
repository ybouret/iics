#include "../visit/simulation.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/string/vfs-utils.hpp"

static inline
Real GetLuaReal( lua_State *L, const char *id )
{
    assert(id);
    lua_settop(L,0);
    lua_getglobal(L, id);
    if( !lua_isnumber(L, -1))
        throw exception("Lua: '%s' is not a number", id);
    return Real(lua_tonumber(L, -1));
}

static inline
bool GetLuaBool( lua_State *L, const char *id )
{
    assert(id);
    lua_settop(L,0);
    lua_getglobal(L, id);
    if( !lua_isboolean(L, -1))
        throw exception("Lua: '%s' is not a boolean", id);
    return lua_toboolean(L,-1) == 1;
}

int main( int argc, char *argv[] )
{
    const char *progname = _vfs::get_base_name(argv[0]);
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
        if( argc <= 1 )
        {
            throw exception("usage: %s config.lua", progname);
        }
        
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
        // Read Lua config
        //----------------------------------------------------------------------
        Lua::State VM;
        lua_State *lua = VM();
        Lua::Config::DoFile(lua, argv[1]);
        
        
        //----------------------------------------------------------------------
        // Setup Simulation
        //----------------------------------------------------------------------
        const Coord  N(30,50);
        const Vertex L(6.0,8.0);
        Simulation   sim(MPI,N,L);
        MPI.PrintfI(stderr, "layout: (%d,%d) -> (%d,%d)\n", int(sim.lower.x),int(sim.lower.y),int(sim.upper.x),int(sim.upper.y));
        MPI.WaitFor(0.1);
      
        SaveGrid( sim.mesh, vformat("g%d.%d.dat",sim.par_size,sim.par_rank));
        
        sim.bubbles.gamma  = 0.01;
        sim.right_pressure = GetLuaReal(lua, "right_pressure");
        sim.right_wall     = GetLuaBool(lua, "right_wall");
        
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