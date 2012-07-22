#include "cell.hpp"

void Cell:: compute_pressure()
{
    const unit_t x1 = lower.x + 1;
    const unit_t xn = upper.x - 1;
    const unit_t y0 = lower.y;
    const unit_t yn = upper.y;
    
    //==========================================================================
    //
    // Fill boundary conditions inside bubble
    //
    //==========================================================================
    for( const Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next )
    {
        for( const GridMarker *g = bubble->markers.head; g; g=g->next )
        {
            P[g->pos.y][g->pos.x] = bubble->pressure;
        }
    }
    
    for( unit_t j=yn;j>=y0;--j )
    {
        P[j][upper.x] = 0.5;
        P[j][lower.x] = 0.5;
    }
    
    for( size_t iter=0; iter < 100; ++iter )
    {
        //======================================================================
        //
        // one red-black iteration
        //
        //======================================================================
        for( unit_t r=0; r <= 1; ++r )
        {
            for( unit_t j=yn;j>=y0;--j )
            {
                const ArrayInt1D &B_j = B[j];
                for( unit_t i=xn-r; i>=x1; i -=2 )
                {
                    if(B_j[i] == 0)
                    {
                        const Real residue = P[j+1][i] + P[j-1][i] + P[j][i+1] + P[j][i-1] - 4 * P[i][j];
                        P[i][j] += 0.25 * residue;
                    }
                }
            }
        }
    }
    
    
}
