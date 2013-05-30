#include "junctions.hpp"

void Junctions:: segment(Array &B) const
{
    assert( grid.is_same_layout_than(B) );
    B.ld(-1);
    
    //==========================================================================
    //
    // Ray Casting from left to right, for each line
    //
    //==========================================================================
    const Real xmin = grid.X()[ grid.lower.x ];
    const Real xmax = grid.X()[ grid.upper.x ];
    assert(xmin<xmax);
    
    for( unit_t y=B.upper.y; y >= B.lower.y; --y)
    {
        
        const Junction::List &JL = Horz(y);
        std::cerr << "\t@y=" << JL.level << " : " << JL.size << std::endl;
        if( JL.size <= 1 )
        {
            //------------------------------------------------------------------
            // 0: no intersection
            // 1: special "tangent" case
            //------------------------------------------------------------------
            continue;
        }
        assert(JL.size>=2);
        
        
    }
    
}
