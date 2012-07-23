#include "yocto/utest/run.hpp"
#include "../segmenter.hpp"


YOCTO_UNIT_TEST_IMPL(seg)
{
    
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
    
    
    { 
        const Layout  lay( Coord(0,0), Coord(20,30) );
        const Region  reg( Vertex(0,pbc.lo/2), Vertex(box.x,pbc.up/2) );
        workspace<Layout, Real, rmesh> W(lay, gs, fs);
        W.mesh.regular_map_to(reg, lay);
        SaveGrid(W.mesh,"grid0.dat");
        Segmenter Seg( W.mesh );
        
        {
            const Vertex dr( Seg.dX[0], Seg.dY[0] );
            bubbles.lambda = dr.norm()/2;
        }
        
        Seg.allocate_segments();
        
        //-- small bubble
        bubbles.create()->map_circle(center, 2);
        bubbles.check_topologies();
        bubbles.first()->save_dat("bubble0.dat");
        
        bubbles.check_geometries_within( reg.vmin.y, reg.vmax.y);

        
        Seg.clear();
        Seg.process_bubble( bubbles.first() );
        Seg.save_junctions("junc0.dat");
        
        //-- bigger bubble
        bubbles.empty();
        
        bubbles.create()->map_circle(center, 3);
        bubbles.check_topologies();
        bubbles.first()->save_dat("bubble1.dat");
        
        bubbles.check_geometries_within( reg.vmin.y, reg.vmax.y);
        
        
        Seg.clear();
        Seg.process_bubble( bubbles.first() );
        Seg.save_junctions("junc1.dat");

    }
    
}
YOCTO_UNIT_TEST_DONE()

