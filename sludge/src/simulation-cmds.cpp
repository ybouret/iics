#include "simulation.hpp"

void Simulation:: perform(const string &cmd)
{
    if( cmd == "save1" )
    {
        if( master )
        {
            bubbles.first()->save_vtk("bubble.vtk");
            bubbles.first()->save_vtk_t("bubble_t.vtk");
            bubbles.first()->save_vtk_n("bubble_n.vtk");
            bubbles.first()->save_dat("bubble.dat");
        }
    }
    
    
    if( cmd == "inter" )
    {
        if( master )
        {
            bubbles.first()->save_dat("bubble.dat");
            bubbles.first()->save_inside("inside.dat", mesh);
            segmenter.save_junctions("junc.dat");
            SaveGrid(mesh, "grid.dat");
        }
        
    }
    
    if( cmd == "raz" )
    {
        //----------------------------------------------------------------------
        // sim initialize P,U,bubbles
        //----------------------------------------------------------------------
        initialize();
        
               
    }
}
