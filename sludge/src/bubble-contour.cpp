#include "bubble.hpp"

void Bubble:: raw_initialize()
{
    
    Tracer *p = root;
    for( size_t i=size;i>0;--i,p=p->next )
    {
        pbc(p->vertex);
        const Tracer *q = p->next; assert(q!=NULL);
        Vertex        pq(p->vertex,q->vertex);
        pbc(pq);
        p->edge = pq;
        p->s2   = p->edge.norm2();
        p->s    = Sqrt( p->s2 );
    }
    compute_area();
    content = pressure * area;
    compute_geometry();
    
}


#if 0
void Bubble:: upgrade_topology()
{
    assert(size>=3);
    assert(root!=NULL);
    
    //--------------------------------------------------------------------------
    // pass 0: pbc on vertices
    //--------------------------------------------------------------------------
    Tracer *p = root;
    
    for( size_t i=size;i>0;--i,p=p->next)
    {
        pbc(p->vertex);
    }
    
   
    
    //--------------------------------------------------------------------------
    // pass : refinement
    //--------------------------------------------------------------------------
    assert(root==p);
    Tracer *q = p->next;
    for( size_t i=0; i < size; ++i )
    {
        
        for(;;)
        {
            //------------------------------------------------------------------
            // compute length to next vertex, with pbc
            //------------------------------------------------------------------
            Vertex pq(p->vertex,q->vertex);
            pbc(pq);
            p->s2 = pq.norm2(); assert(p->s2>0);
            p->s  = Sqrt( p->s2 );
            
            //------------------------------------------------------------------
            // do we refine ?
            //------------------------------------------------------------------
            if( p->s > lambda )
            {
                Tracer *I   = cache.provide();
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
            p->edge= pq;
            assert(p->s2>0);
            assert(p->s>0);
            break;
        }
        
        //----------------------------------------------------------------------
        // next edge
        //----------------------------------------------------------------------
        p = q;
        q = q->next;
    }
    

}
#endif
