#include "junctions.hpp"
#include "yocto/exception.hpp"

void Junctions:: segment(Array &B) const
{
    assert( grid.is_same_layout_than(B) );
    B.ld(-1);
    
    //==========================================================================
    //
    // Ray Casting from left to right, for each line
    //
    //==========================================================================
    const Array1D &X    = grid.X();
    const Real     xmin = X[ grid.lower.x ];
    const Real     xmax = X[ grid.upper.x ];
    assert(xmin<xmax);
    
    for( unit_t y=B.upper.y; y >= B.lower.y; --y)
    {
        
        const Junction::List &JL = Horz(y);
        
        if(JL.size<=1)
            continue;
        
        
        bool  in_bubble= true;
        const Junction *J = JL.head; assert(J);
        const Junction *K = J->next;
        while(K)
        {
            if(in_bubble)
            {
                if(J->owner!=K->owner)
                    throw exception("Bubble in Bubble!");
                if(J->inside||K->inside)
                {
                    const unit_t ini = J->inside ? J->upper : B.lower.x;
                    const unit_t end = K->inside ? K->lower : B.upper.x;
                    if(end>=ini)
                    {
                        const Real u = J->owner->UID;
                        for(unit_t i=ini;i<=end;++i)
                            B[y][i] = u;
                    }
                }
            }
            in_bubble = !in_bubble;
            J=K;
            K=K->next;
        }
        
    }
    
}


void Junctions:: save_inside_of( const Array &B, const string &fn ) const
{
    assert(B.is_same_layout_than(grid) );
    ios::ocstream fp(fn,false);
    for(unit_t j=B.lower.y;j<=B.upper.y;++j)
    {
        for(unit_t i=B.lower.x;i<=B.upper.x;++i)
        {
            if(B[j][i]>=0)
                fp("%g %g\n", grid.X()[i], grid.Y()[j]);
        }
    }
}
