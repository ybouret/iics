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
        
        {
            ios::ocstream fp("grid.dat",false);
            for( unit_t j=sim.lower.y; j<=sim.upper.y;++j)
            {
                if( (j&1) )
                {
                    for( unit_t i=sim.lower.x; i <= sim.upper.x; ++i )
                    {
                        fp("%g %g\n", sim.X[i], sim.Y[j] );
                    }
                }
                else
                {
                    for( unit_t i=sim.upper.x; i >= sim.lower.x; --i )
                    {
                        fp("%g %g\n", sim.X[i], sim.Y[j] );
                    }
                }
            }
            
            for( unit_t i=sim.lower.x; i <= sim.upper.x; ++i )
            {
                if( (i&1) )
                {
                    for( unit_t j=sim.lower.y; j<=sim.upper.y;++j)
                    {
                        fp("%g %g\n", sim.X[i], sim.Y[j] );
                    }
                }
                else {
                    for( unit_t j=sim.upper.y; j>=sim.lower.y;--j)
                    {
                        fp("%g %g\n", sim.X[i], sim.Y[j] );
                    }
                }
            }
        }
        
        //----------------------------------------------------------------------
        // sim initialize P,U,bubbles
        //----------------------------------------------------------------------
        sim.initialize();
        
        //----------------------------------------------------------------------
        // Initialize bubbles
        //----------------------------------------------------------------------
        sim.check_and_dispatch_bubbles();
        
        
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