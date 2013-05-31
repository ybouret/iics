#include "workspace.hpp"

void Workspace:: reset_pressure()
{
    P.ldz();
    pressurize();
}

void Workspace:: pressurize()
{
    //--------------------------------------------------------------------------
    // Fill the bubbles location with bubble pressures
    //--------------------------------------------------------------------------
    bpres.free();
    bpres.reserve(bubbles.size);
    for( const Bubble *b = bubbles.head;b;b=b->next)
    {
        bpres.push_back( b->pressure );
    }
    
    for( unit_t j=outline.lower.y; j <= outline.upper.y; ++j)
    {
        for(unit_t i=lower.x;i<=upper.x;++i)
        {
            const Real which = B[j][i];
            if( which >= 0 )
            {
                assert(which<bubbles.size);
                assert(which<bpres.size());
                P[j][i] = bpres[ which+1 ];
            }
        }
    }
    
    
}