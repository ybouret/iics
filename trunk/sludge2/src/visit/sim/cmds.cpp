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
        bubbles.first()->save_vtk_gn( "bgn.vtk" );
    }
}

