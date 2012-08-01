#include "yocto/utest/run.hpp"
#include "../bubble.hpp"
#include "../rescaler.hpp"

static inline void save_all( const Bubble &bubble, const string &pfx )
{
    bubble.save_dat(   pfx + ".dat" );
    bubble.save_vtk(   pfx + ".vtk");
    bubble.save_vtk_t( pfx + "_t.vtk");
    bubble.save_vtk_n( pfx + "_n.vtk");
}

static inline
void expand( Bubble &bubble , const Vertex &center)
{
    Tracer *p = bubble.root;
    for( size_t i=bubble.size;i>0;--i,p=p->next)
    {
        const Real    theta = p->vertex.angle();
        const Vertex  v( 1.5 * Cos(theta), 2.0 * Cos( theta ) );
        const Vertex  r = p->vertex - center;
        p->vertex = center + v.norm() * r;
    }
}

YOCTO_UNIT_TEST_IMPL(shapes)
{
    
    double lambda = 1;
    Vertex box(100,100);
    PBC    pbc(box.y);
    Tracer::Cache tcache;
    Spot::Cache   scache;
    Marker::Cache mcache;
    Rescaler      rescaler;
    Vertex center( box.x/2, 0.0 );
    
    Bubble bubble(lambda,pbc,tcache,scache,mcache);
#if 1
    bubble.map_circle( center, 2.0 + Alea() );
    std::cerr << "curvature: " << bubble.root->curvature << std::endl;
    save_all( bubble, "circle" );
    rescaler.upgrade(bubble);
    bubble.compute_geometry();
    bubble.save_dat("circle1.dat");
    
    expand(bubble,center);
    std::cerr << "Rescaling circle" << std::endl;
    rescaler.upgrade(bubble);
    bubble.save_dat("circle2.dat");
#endif
    
    bubble.map_peanut( center, 2.0+Alea(), 0.96);
    save_all( bubble, "peanut");
    rescaler.upgrade(bubble);
    bubble.save_dat("peanut1.dat");
    
    expand(bubble,center);
    rescaler.upgrade(bubble);
    bubble.save_dat("peanut2.dat");

}
YOCTO_UNIT_TEST_DONE()
