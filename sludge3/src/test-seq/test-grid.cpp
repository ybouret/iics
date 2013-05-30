#include "yocto/utest/run.hpp"
#include "yocto/code/rand.hpp"
#include "junctions.hpp"
#include "shape.hpp"

static inline
void perform_locate( const Grid &grid, const Bubble &bubble )
{
    const Tracer *tr = bubble.root;
    for(size_t i=bubble.size;i>0;--i,tr=tr->next)
    {
        Coord C;
        const int ans = __Grid::Locate(grid, tr->pos, C);
        if( SLUDGE_INSIDE == ans)
        {
        std::cerr << tr->pos << " : " << ans << " @" << C << std::endl;
        }
        else
            std::cerr << tr->pos << " : outside=" << ans << std::endl;
    }
}

#include "yocto/string/conv.hpp"

YOCTO_UNIT_TEST_IMPL(grid)
{
    
    Real radius = 5.1;
    if( argc > 1 )
    {
        radius = strconv::to<Real>( argv[1], "radius" );
    }
    
    std::cerr << "*** Building grid" << std::endl;
    array_db       adb;
    const Layout   L( Coord(1,1), Coord(10,20) );
    Grid           grid(L,adb);
    const Region2D R( Vertex(-5,-6), Vertex(5,6) );
    
    Junctions junctions(grid);
    grid.regular_map_to(R,L);
    
    __Grid::SaveDat(grid, "grid.dat" );
    
    std::cerr << "*** Building bubble" << std::endl;
    Real lam = __Grid::ComputeLambda(grid);
    Real gam = 1;
    std::cerr << "lambda=" << lam << std::endl;
    
    Bubble bubble(lam,gam,1);
    Shape::Circle(&bubble, Vertex(0,0), radius);
    bubble.regularize();
    bubble.save_all("b1");
   
    std::cerr << "*** Testing bubble locations" << std::endl;
    perform_locate(grid, bubble);
    
    for(unit_t i=grid.lower.x; i<= grid.upper.x; ++i)
    {
        const Junction::List &J = junctions.Vert(i);
        assert( Junction::Vert == J.type);
        std::cerr << "Vert(" << i << ")@" << J.level << std::endl;
    }
    
    for(unit_t j=grid.lower.y; j<= grid.upper.y; ++j)
    {
        const Junction::List &J = junctions.Horz(j);
        assert( Junction::Horz == J.type);
        std::cerr << "Horz(" << j << ")@" << J.level << std::endl;
    }
    
    std::cerr << "*** Segmenting bubble" << std::endl;
    junctions.clear();
    bubble.compute_curvatures();
    junctions.inter(bubble);
    junctions.sort();
    junctions.save_dat( "j1.dat" );
    junctions.save_t( "j1_t.dat");
    junctions.save_n( "j1_n.dat");

    Shape::Blob(&bubble, Vertex(0,0), radius, 0.7, 0.6);
    bubble.regularize();
    bubble.save_all("b2");
    junctions.clear();
    junctions.inter(bubble);
    junctions.sort();
    junctions.save_dat( "j2.dat" );
    junctions.save_t("j2_t.dat");
    junctions.save_n("j2_n.dat");

}
YOCTO_UNIT_TEST_DONE()

