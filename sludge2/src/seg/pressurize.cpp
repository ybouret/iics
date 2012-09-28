#include "../segmenter.hpp"

void Segmenter:: pressurize(Array &P) const
{
    
    for( const Marker *m = markers.head;m;m=m->next)
    {
        const Bubble *bubble = m->bubble; assert(bubble);
        P[m->inside.y][m->inside.x] = bubble->pressure;
    }
    
}


void Segmenter:: build_effective_pressure( const Array &B, VertexArray &Penter, VertexArray &Pleave )
{
    Penter.ldz();
    Pleave.ldz();

#if 1
    for( const Marker *m = markers.head;m;m=m->next)
    {
        const Bubble *bubble = m->bubble; assert(bubble);
        Vertex &Pe = Penter[m->inside.y][m->inside.x];
        Pe.x = Pe.y = bubble->pressure;
        Vertex &Pl = Pleave[m->inside.y][m->inside.x];
        Pl.x = Pl.y = bubble->pressure;
    }
#endif
    
    const double gamma = 1;
    //--------------------------------------------------------------------------
    // horizontal setting
    //--------------------------------------------------------------------------
    for( unit_t j=Y.lower;j<=Y.upper;++j)
    {
        const Segment &seg     = Horz(j);
        const Junction *J      = seg.head;
        bool            inside = false;
        while(J)
        {
            size_t count = 1;
            const Junction *K  = J; //!< first junction @J->klo
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
                const Bubble *bubble = J->bubble; assert(bubble);
                if(inside)
                {
                    //! leaving the bubble: take the last curvature
                    assert( B[j][J->klo] >0 );
                    Pleave[j][J->klo].x = bubble->pressure + gamma * J->curvature;
                    /*
                     
                     while(curr<=J->klo)
                     {
                     //B_j[curr]   = bubble->id;
                     //P_j[curr].x = bubble->pressure;
                     Marker *m   = markers.append();
                     m->inside.x = curr;
                     m->inside.y = j;
                     m->bubble   = bubble;
                     ++curr;
                     }
                     */
                }
                else
                {
                    //! entering the bubble: take the first curvature
                    assert(B[j][K->klo]<=0);
                    Penter[j][K->khi].x = bubble->pressure + gamma * K->curvature;
                }
                inside = !inside;
            }
            //curr = J->khi;
            J    = J->next;
        }
        
    }
    
    
}
