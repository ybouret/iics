#include "../cell.hpp"




int main( int argc, char *argv[] )
{
    
    try 
    {
        mpi &MPI = mpi::init( &argc, &argv);
        
        Cell cell(10,20,5.0,6.0,MPI);
        MPI.Printf(stderr, "rank %d> SubRegion y: %g -> %g\n", MPI.CommWorldRank, cell.SubRegion.vmin.y, cell.SubRegion.vmax.y );
        
        
        if( MPI.IsMaster )
        {
            Bubble *b = cell.bubbles.create();
            b->map_circle(V2D(0,0), 1.2);
        }
        
        cell.bubbles.dispatch_all(MPI);
        
        
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "Unhandled exception" << std::endl;
    }
    return 1;
}
