#include "simulation.hpp"

static 
void save_inside( const array<V2D> &pts )
{
    ios::ocstream fp( "inside.dat", false );
    for( size_t i=pts.size(); i>0; --i )
    {
        fp("%g %g\n",pts[i].x, pts[i].y);
    }
}

int main( int argc, char *argv[] )
{
    
    try 
    {
        AleaInit();
        
        //----------------------------------------------------------------------
        //info for VisIt
        //----------------------------------------------------------------------
        const string sim_name    = "Laponite";
        const string sim_comment = "Laponite/Hele-Shaw";
        const string sim_path    = ".";
        
        //----------------------------------------------------------------------
        // initialize MPI
        //----------------------------------------------------------------------
        mpi &MPI = mpi::init( &argc, &argv );
        const int sim_rank = MPI.CommWorldRank;
        //const int sim_size = MPI.CommWorldSize;
        
        const string       trace_name = "trace.dat";
        VisIt:: TraceFile  trace_file(sim_rank,trace_name );
        
        //----------------------------------------------------------------------
        // parallel VisIt
        //----------------------------------------------------------------------
        VisIt:: SetupParallel( MPI, sim_name, sim_comment, sim_path, NULL);
        
        Simulation sim(15,18,5.0,6.0,MPI);
        
        //----------------------------------------------------------------------
        // sim initialize P,U,bubbles
        //----------------------------------------------------------------------
        sim.initialize();
        
        //----------------------------------------------------------------------
        // Initialize bubbles
        //----------------------------------------------------------------------
        sim.check_and_dispatch_bubbles();
        vector<V2D> pts;
        
        sim.collect_inside(pts);
        sim.bubbles.first()->save_dat("bubble.dat");
        save_inside(pts);
        
        
        
        return 0;
        
        //----------------------------------------------------------------------
        // Initialize ghosts
        //----------------------------------------------------------------------
        sim.init_exchange();
        sim.wait_exchange();
        
        
        //----------------------------------------------------------------------
        // call VisIt
        //----------------------------------------------------------------------
        VisIt:: MainLoop(sim);
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "Unhandled exception" << std::endl;
    }
    return 1;
}