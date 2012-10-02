#include "../segmenter.hpp"

#if 0
void Segmenter:: pressurize(Array &P) const
{
    
    for( const Marker *m = markers.head;m;m=m->next)
    {
        const Bubble *bubble = m->bubble; assert(bubble);
        P[m->inside.y][m->inside.x] = bubble->pressure;
    }
    
}
#endif

void Segmenter:: build_effective_pressure( const Array &B, Array &P, VertexArray &Penter, VertexArray &Pleave )
{
    Penter.ldz();
    Pleave.ldz();
    
    for( const Marker *m = markers.head;m;m=m->next)
    {

        const Bubble *bubble = m->bubble; assert(bubble);
        P[m->inside.y][m->inside.x] = bubble->pressure;
#if 0
        Vertex &Pe = Penter[m->inside.y][m->inside.x];
        Pe.x = Pe.y = bubble->pressure;
        Vertex &Pl = Pleave[m->inside.y][m->inside.x];
        Pl.x = Pl.y = bubble->pressure;
#endif
    }
    
    const Real gamma = 0.0;
    
    //--------------------------------------------------------------------------
    //
    // horizontal effective pressure
    //
    //--------------------------------------------------------------------------
    for( unit_t j=Y.lower;j<=Y.upper;++j)
    {
        const Segment &seg     = Horz(j);
        const Junction *J      = seg.head;
        
        while(J)
        {
            size_t          count = 1;
            const Junction *K     = J; //!< first junction @J->klo
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
                //--------------------------------------------------------------
                // we crossed a bubble: check K
                //--------------------------------------------------------------
                if( B[j][K->klo] > 0 )
                {
                    assert(B[j][J->khi]<=0);
                    //----------------------------------------------------------
                    // we leave a bubble
                    //----------------------------------------------------------
                    Pleave[j][J->klo].x = J->bubble->pressure + gamma * J->curvature;
                }
                else
                {
                    assert(B[j][J->khi]>0);
                    //----------------------------------------------------------
                    // we enter a bubble
                    //----------------------------------------------------------
                    Penter[j][K->khi].x = K->bubble->pressure + gamma * K->curvature;
                }
            }
            J    = J->next;
        }
    }
    
    ios::ocstream fp("jvert.dat",false);
    //--------------------------------------------------------------------------
    //
    // vertical effective pressure
    //
    //--------------------------------------------------------------------------
    for( unit_t i=X.lower;i<=X.upper;++i)
    {
        const Segment &seg     = Vert(i);
        const Junction *J      = seg.head;
        fp("@i=%d\n",i);
        while(J)
        {
            size_t          count = 1;
            const Junction *K     = J; //!< first junction @J->klo
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
                //--------------------------------------------------------------
                // we crossed a bubble: check K
                //--------------------------------------------------------------
                if( B[K->klo][i] > 0 )
                {
                    assert(B[J->klo][i]>0);
                    assert(B[J->khi][i]<=0);
                    //----------------------------------------------------------
                    // we leave a bubble
                    //----------------------------------------------------------
                    Pleave[J->klo][i].y = J->bubble->pressure + gamma * J->curvature;
                    fp("Pleave[%d][%d].y=%g\n", J->klo, i, Pleave[J->klo][i].y );
                }
                else
                {
                    assert(B[K->klo][i]<=0);
                    assert(B[K->khi][i]>0);
                    //----------------------------------------------------------
                    // we enter a bubble
                    //----------------------------------------------------------
                    Penter[K->khi][i].y = K->bubble->pressure + gamma * K->curvature;
                    fp("Penter[%d][%d].y=%g\n", K->khi, i, Penter[K->khi][i].y );

                }

            }
            J=J->next;
        }
    }
    
    
    
}
