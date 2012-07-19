#include "bubble.hpp"

void Bubble:: update_contour()
{
    assert(size>=3);
    assert(root!=NULL);
    
    //--------------------------------------------------------------------------
    // pass 0: pbc on vertices
    //--------------------------------------------------------------------------
    Point *p = root;
    
    for( size_t i=size;i>0;--i,p=p->next)
        pbc(p->vertex);
    
    p        = root;
    Point *q = p->next;
    
    //--------------------------------------------------------------------------
    // pass : refinement
    //--------------------------------------------------------------------------
    for( size_t i=0; i < size; ++i )
    {
        
        for(;;)
        {
            //------------------------------------------------------------------
            // compute length to next vertex, with pbc
            //------------------------------------------------------------------
            V2D pq(p->vertex,q->vertex);
            pbc(pq);
            p->s_next = pq.norm(); assert(p->s_next>0);
            
            //------------------------------------------------------------------
            // do we refine ?
            //------------------------------------------------------------------
            if( p->s_next > lambda )
            {
                Point *I = create();
                I->vertex.x = p->vertex.x + 0.5 * pq.x;
                I->vertex.y = p->vertex.y + 0.5 * pq.y;
                pbc(I->vertex);
                insert_after(p, I);
                assert(p->next==I);
                assert(I->prev==p);
                q = I;
                continue;
            }
            
            //------------------------------------------------------------------
            // ok, keep that in mind...
            //------------------------------------------------------------------
            p->r_next = pq;
            assert(p->s_next>0);
            break;
        }
        
        //----------------------------------------------------------------------
        // next edge
        //----------------------------------------------------------------------
        p = q;
        q = q->next;
    }
    
    p = root;
    for( size_t i=0; i < size; ++i, p=p->next)
    {
        assert(p->s_next>0);
    }
    
}
