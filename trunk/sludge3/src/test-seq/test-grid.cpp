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

YOCTO_UNIT_TEST_IMPL(grid)
{
    array_db       adb;
    const Layout   L( Coord(1,1), Coord(10,20) );
    Grid           grid(L,adb);
    const Region2D R( Vertex(-5,-6), Vertex(5,6) );
    
    Junctions junctions(grid);
    grid.regular_map_to(R,L);
    
    __Grid::SaveDat(grid, "grid.dat" );
    
    Real lam = __Grid::ComputeLambda(grid);
    std::cerr << "lambda=" << lam << std::endl;
    
    Bubble bubble(lam);
    Shape::Circle(&bubble, Vertex(0,0), 5.1);
    bubble.auto_contour();
    bubble.save_dat("bubble.dat");
   
    perform_locate(grid, bubble);
    
    for(unit_t i=junctions.lower.x; i<= junctions.upper.x; ++i)
    {
        const Junction::List &J = junctions.Vert(i);
        assert( Junction::Vert == J.type);
        std::cerr << "Vert(" << i << ")@" << J.level << std::endl;
    }
    
    for(unit_t j=junctions.lower.y; j<= junctions.upper.y; ++j)
    {
        const Junction::List &J = junctions.Horz(j);
        assert( Junction::Horz == J.type);
        std::cerr << "Horz(" << j << ")@" << J.level << std::endl;
    }
    
}
YOCTO_UNIT_TEST_DONE()

