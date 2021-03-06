#include "../segmenter.hpp"

void Segmenter:: build_bubbles_field( Array &B )
{
    B.ldz();
    assert( 0 == markers.size);
    
    //fprintf(stderr, "[Bubble Field]\n");
    //--------------------------------------------------------------------------
    // scan/line
    //--------------------------------------------------------------------------
    for(unit_t j=Y.lower;j<=Y.upper;++j)
    {
        const Segment &seg     = Horz(j);
        const Junction *J      = seg.head;
        unit_t          curr   = X.lower;
        bool            inside = false;
        Array1D        &B_j    = B[j];
        //if( seg.size ) fprintf( stderr, "@j=%ld: Y=%g\n", j, Y[j]);
        
        while(J)
        {
            //------------------------------------------------------------------
            // we met a junction
            // count how many junctions between two vertices
            //------------------------------------------------------------------
            size_t count = 1;
            while(J->next &&
                  (J->next->klo == J->klo) )
            {
                assert(J->bubble);
                assert(J->next->bubble);
                if(J->bubble != J->next->bubble)
                    throw exception("Invalid overlapping bubbles in Segmenter::build");
                J=J->next;
                ++count;
            }
            //fprintf( stderr, "\tstatus=%s count=%lu @klo=%ld,khi=%ld\n", (inside ? "inside " : "outside"), count, J->klo, J->khi);
            
            //------------------------------------------------------------------
            // update current values
            //------------------------------------------------------------------
            if(inside)
            {
                while(curr<=J->klo)
                {
                    const Bubble *bubble = J->bubble; assert(bubble);
                    B_j[curr]   = bubble->id;
                    Marker *m   = markers.append();
                    m->inside.x = curr;
                    m->inside.y = j;
                    m->bubble   = bubble;
                    ++curr;
                }
            }
            
            //------------------------------------------------------------------
            // update current status
            //------------------------------------------------------------------
            if(count&1)
            {
                inside = !inside;
            }
            
            curr = J->khi;
            J    = J->next;
        }
    }
    
}


