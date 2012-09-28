#include "../segmenter.hpp"

void Segmenter:: pressurize(Array &P) const
{

    for( const Marker *m = markers.head;m;m=m->next)
    {
        const Bubble *bubble = m->bubble; assert(bubble);
        P[m->inside.y][m->inside.x] = bubble->pressure;
    }

}


void Segmenter:: build_effective_pressure( VertexArray &Peff )
{
    Peff.ldz();
    //--------------------------------------------------------------------------
    // horizontal setting
    //--------------------------------------------------------------------------
    for( unit_t j=Y.lower;j<=Y.upper;++j)
    {
        const Segment &seg     = Horz(j);
        const Junction *J      = seg.head;
        unit_t          curr   = X.lower;
        bool            inside = false;
        VertexArray1D  &P_j    = Peff[j];
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
                    //! leaving the bubble
                    while(curr<=J->klo)
                    {
                        const Bubble *bubble = J->bubble; assert(bubble);
                        //B_j[curr]   = bubble->id;
                        //P_j[curr].x = bubble->pressure;
                        Marker *m   = markers.append();
                        m->inside.x = curr;
                        m->inside.y = j;
                        m->bubble   = bubble;
                        ++curr;
                    }
                }
                else
                {
                    //! entering the bubble
                }
                inside = !inside;
            }
            curr = J->khi;
            J    = J->next;
        }

    }
    
    
}
