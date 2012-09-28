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
    
#if 0
    for( const Marker *m = markers.head;m;m=m->next)
    {
        const Bubble *bubble = m->bubble; assert(bubble);
        Vertex &Pe = Penter[m->inside.y][m->inside.x];
        Pe.x = Pe.y = bubble->pressure;
        Vertex &Pl = Pleave[m->inside.y][m->inside.x];
        Pl.x = Pl.y = bubble->pressure;
    }
#endif
    
    const double gamma = 0.1;
    //--------------------------------------------------------------------------
    //
    // horizontal effective pressure
    //
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
                    throw exception("Invalid overlapping bubbles in Segmenter::build_effective_pressure/horz");
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
                }
                else
                {
                    //! entering the bubble: take the first curvature
                    assert(B[j][K->klo]<=0);
                    Penter[j][K->khi].x = bubble->pressure + gamma * K->curvature;
                }
                inside = !inside;
            }
            J    = J->next;
        }
    }
    
    return;
    
    //--------------------------------------------------------------------------
    //
    // vertical effective pressure
    //
    //--------------------------------------------------------------------------
    for( unit_t i=X.lower;i<=X.upper;++i)
    {
        const Segment &seg     = Vert(i);
        const Junction *J      = seg.head;
        bool            inside = B[Y.lower][i] > 0;
        while(J)
        {
            std::cerr << "J->klo=" << J->klo << "/Y.lower=" << Y.lower << std::endl;
            size_t count = 1;
            const Junction *K  = J; //!< first junction @J->klo
            while(J->next && J->next->klo == J->klo )
            {
                assert(J->bubble);
                assert(J->next->bubble);
                if(J->bubble != J->next->bubble)
                    throw exception("Invalid overlapping bubbles in Segmenter::build_effective_pressure/vert");
                J=J->next;
                ++count;
            }
            if(count&1)
            {
                const Bubble *bubble = J->bubble; assert(bubble);
                if(inside)
                {
                    //! leaving the bubble: take the last curvature
                    assert( B[J->klo][i]>0);
                    Pleave[J->klo][i].y = bubble->pressure + gamma * J->curvature;
                }
                else
                {
                    //! entering the bubble: take the first curvature
                    assert(B[K->klo][i]<=0);
                    Penter[K->khi][i].y = bubble->pressure + gamma * K->curvature;
                }
                inside = !inside;
            }
            J    = J->next;
        }
    }

    
    
}
