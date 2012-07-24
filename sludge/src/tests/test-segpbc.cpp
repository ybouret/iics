#include "yocto/utest/run.hpp"
#include "../segmenter.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/swamp/vtk-writer.hpp"

static inline void __process(Bubbles              &bubbles, 
                             Segmenter            &Seg, 
                             int level,
                             const Region         &reg,
                             vtk_writer           &vtk,
                             const WorkspaceBase  &W,
                             Array                &B)
{
    bubbles.check_topologies();
    bubbles.first()->save_dat( vformat("bubble%d.dat",level) );
    bubbles.check_geometries_within( reg.vmin.y, reg.vmax.y);
    
    
    Seg.process_bubbles( bubbles );
    Seg.horizontal_pbc();
    Seg.save_junctions( vformat("junc%d.dat",level) );
    Seg.assign_markers();
    
    bubbles.fill( B );
    
    vector<string> var;
    var.push_back( "B" );
    vtk.save( vformat("b%d.vtk",level), vformat("b%d",level), W, var, W.__layout());
    bubbles.first()->save_inside( vformat("inside%d.dat",level), W.mesh);
}


YOCTO_UNIT_TEST_IMPL(segpbc)
{
    AleaInit();
    //==========================================================================
    //
    // prepare the grid
    // 
    //==========================================================================
    const Vertex  box(10,10);
    const Vertex  center( box.x/2, 0);
    Bubbles       bubbles( box );
    const PBC    &pbc = bubbles.pbc;
    
    GhostsSetup   gs;
    FieldsSetup   fs;
    
    Y_SWAMP_DECL_AUX(fs, "B", Array);
    
    const Layout  lay( Coord(0,0), Coord(20,30) );
    const Region  reg( Vertex(0,pbc.lo), Vertex(box.x,pbc.up) );
    workspace<Layout, Real, rmesh> W(lay, gs, fs);
    W.mesh.regular_map_to(reg, lay);
    SaveGrid(W.mesh,"grid0.dat");
    Segmenter Seg( W.mesh );
    Array &B = W["B"].as<Array>();
    vtk_writer vtk;
    
    
    bubbles.lambda = min_of( Seg.dX[0], Seg.dY[0])/2;
    
    Seg.allocate_segments();
    
    //-- make a peanut
    bubbles.empty();
    
    bubbles.create()->map_peanut(center+Vertex(0,-4), 3.5, 0.95);
    
    for( int i=0; i < 250; ++i )
    {
        __process(bubbles, Seg, i, reg, vtk, W, B);
        bubbles.first()->translate( Vertex(0,-0.1+0.02*Alea()) );
    }
    
    
}
YOCTO_UNIT_TEST_DONE()
