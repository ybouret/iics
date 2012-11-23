#include "../simulation.hpp"


void Simulation:: perform(const string &cmd, const array<string> &args)
{
    if( cmd == "raz" )
    {
        initialize();
        return;
    }
    
    if( cmd == "geo" )
    {
        segmenter.save_vtk_gn( "jgn.vtk" );
        segmenter.save_vtk_gt( "jgt.vtk" );
        bubbles.first()->save_vtk_gn( "bgn.vtk" );
        bubbles.first()->save_vtk_gt( "bgt.vtk" );
        bubbles.first()->save_vtk( "b.vtk" );
        return;
    }
    
    std::cerr << "\t\t<<Unknown Command>>" << std::endl;
}

