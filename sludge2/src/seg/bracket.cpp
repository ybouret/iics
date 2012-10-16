#include "../segmenter.hpp"



void Segmenter:: find_bracketing_junctions(ConstJunctionPtr &jprev,
                                           ConstJunctionPtr &jnext,
                                           const Spot       *spot) const
{
    assert(spot);
    assert(0==jnext);
    assert(0==jprev);
    const Tracer *tracer = spot->handle; assert(tracer);
    //--------------------------------------------------------------------------
    // find next
    //--------------------------------------------------------------------------
    {
        const Tracer *p  = tracer;
        if( p->jnext )
        {
            jnext = p->jnext;
        }
        else
        {
            for(p=p->next;p!=tracer;p=p->next)
            {
                if(p->jprev)
                {
                    jnext = p->jprev;
                    break;
                }
                
                if(p->jnext)
                {
                    jnext = p->jnext;
                    break;
                }
            }
        }
        
    }
    
    if(!jnext)
    {
        throw exception("can't find jnext for tracer@(%g,%g)!",tracer->vertex.x,tracer->vertex.y);
    }
    
    //--------------------------------------------------------------------------
    // find prev
    //--------------------------------------------------------------------------
    {
        const Tracer *p = tracer;
        if( p->jprev )
        {
            jprev = p->jprev;
        }
        else
        {
            for(p=p->prev;p!=tracer;p=p->prev)
            {
                if(p->jnext)
                {
                    jprev = p->jnext;
                    break;
                }
                
                if(p->jprev)
                {
                    jprev = p->jprev;
                    break;
                }
            }
        }
    }
    if(!jprev)
    {
        throw exception("can't find jprev for tracer@(%g,%g)!",tracer->vertex.x,tracer->vertex.y);

    }
    

    if(jprev==jnext)
    {
        throw exception("same next/prev junction for tracer@(%g,%g)!",tracer->vertex.x,tracer->vertex.y);
    }
    
}
