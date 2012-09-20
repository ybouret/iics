#include "../segmenter.hpp"

void Segmenter:: build( Array &B ) const
{
    B.ldz();
    
       
    //--------------------------------------------------------------------------
    // scan/line
    //--------------------------------------------------------------------------
    for(unit_t j=Y.lower;j<=Y.upper;++j)
    {
        const Segment &seg   = Horz(j);
        const Junction *J    = seg.head;
        unit_t          curr   = X.lower;
        bool            inside = false;
        Array1D        &B_j    = B[j];
        while(J)
        {
            size_t count = 1;
            while(J->next && J->next->klo == J->klo )
            {
                assert(J->bubble);
                assert(J->next->bubble);
                if(J->bubble != J->next->bubble)
                    throw exception("Invalid overlapping bubbles in Segmenter::build");
                J=J->next;
                ++count;
            }
            if(count&1)
            {
                if(inside)
                {
                    while(curr<=J->klo)
                    {
                        B_j[curr++] = J->bubble->id;
                    }
                }
                inside = !inside;
            }
            curr = J->khi;
            J    = J->next;
        }
    }
    
}
