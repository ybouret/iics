#include "../segmenter.hpp"


void Segmenter:: build_effective_pressure(const Array  &B,
                                          Array        &P,
                                          VertexArray  &Penter,
                                          VertexArray  &Pleave,
                                          const Real    gamma)
{
    Penter.ldz();
    Pleave.ldz();
    
    for( const Marker *m = markers.head;m;m=m->next)
    {
        
        const Bubble *bubble = m->bubble; assert(bubble);
        P[m->inside.y][m->inside.x] = bubble->pressure;
    }
    
    save( "j.dat" );
    
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
            size_t          count = 1; //!< one junction => one cross
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
            assert(J->klo==K->klo);
            assert(J->khi==K->khi);
            if(count&1)
            {
                //--------------------------------------------------------------
                // we crossed a bubble: check K
                //--------------------------------------------------------------
                if( B[j][K->klo] > 0 )
                {
                    //----------------------------------------------------------
                    // we leave a bubble: take the last junction
                    //----------------------------------------------------------
                    assert(B[j][J->khi]<=0);
                    Pleave[j][J->klo].x = J->pressure;
                }
                else
                {
                    //----------------------------------------------------------
                    // we enter a bubble: take the first junction
                    //----------------------------------------------------------
                    if( B[j][K->khi] <= 0 )
                    {
                        fprintf( stderr, "error @j=%ld, K->khi=%ld\n and B[klo]=%g @(%g,%g)\n", j, K->khi, B[j][K->klo], X[K->klo],Y[j] );
                    }
                    assert(B[j][K->khi]>0);
                    Penter[j][K->khi].x = K->pressure;
                }
            }
            J    = J->next;
        }
    }
    
    //--------------------------------------------------------------------------
    //
    // vertical effective pressure
    // restriction because of PBC
    //
    //--------------------------------------------------------------------------
    const unit_t ymin=Y.lower + 1;
    const unit_t ymax=Y.upper - 1;
    
    for( unit_t i=X.lower;i<=X.upper;++i)
    {
        const Segment &seg     = Vert(i);
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
                    throw exception("Invalid overlapping bubbles in Segmenter::build_effective_pressure/vert");
                J=J->next;
                ++count;
            }
            assert(J->klo==K->klo);
            assert(J->khi==K->khi);
            
            // restriction
            if((count&1) &&
               (J->klo >= ymin) &&
               (J->khi <= ymax)
               )
            {
                //--------------------------------------------------------------
                // we crossed a bubble: check K
                //--------------------------------------------------------------
                if( B[K->klo][i] > 0 )
                {
                    //----------------------------------------------------------
                    // we leave a bubble
                    //----------------------------------------------------------
                    assert(B[J->klo][i]>0);
                    assert(B[J->khi][i]<=0);
                    Pleave[J->klo][i].y = J->pressure;
                }
                else
                {
                    //----------------------------------------------------------
                    // we enter a bubble
                    //----------------------------------------------------------
#if 0
                    if( !(B[K->khi][i]>0) )
                    {
                        fprintf( stderr, "Error khi=%ld,i=%ld (x=%g,y=%g)\n",K->khi,i,X[i],Y[K->khi]);
                    }
#endif
                    assert(B[K->klo][i]<=0);
                    assert(B[K->khi][i]>0);
                    
                    Penter[K->khi][i].y = K->pressure;
                }
            }
            
            J=J->next;
        }
    }
}
