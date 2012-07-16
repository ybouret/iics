#include "../cell.hpp"




int main( int argc, char *argv[] )
{
    
    try 
    {
        mpi &MPI = mpi::init( &argc, &argv);
        
        Parameters param(10,20,15.0,30.0,MPI);
        MPI.Printf(stderr, "rank %d> SubRegion y: %g -> %g\n", MPI.CommWorldRank, param.SubRegion.vmin.y,param.SubRegion.vmax.y );
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "Unhandled exception" << std::endl;
    }
    return 1;
}
