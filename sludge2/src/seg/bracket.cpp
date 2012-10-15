#include "../segmenter.hpp"



void Segmenter:: find_bracketing_junctions(const Spot *spot) const
{
    assert(spot);
    
    const Tracer *tracer = spot->handle; assert(tracer);
    //fprintf( stderr, "\ttracer @(%g,%g)\n", tracer->vertex.x, tracer->vertex.y);
    //--------------------------------------------------------------------------
    // find next
    //--------------------------------------------------------------------------
    const Junction *jnext = 0;
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
        fprintf( stderr, "can't find jnext!\n");
    }
    //fprintf( stderr, "\t\t: jnext@(%g,%g)\n", jnext->vertex.x, jnext->vertex.y);
    
    //--------------------------------------------------------------------------
    // find prev
    //--------------------------------------------------------------------------
    const Junction *jprev = 0;
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
        fprintf( stderr, "can't find jprev!");
    }
    //fprintf( stderr, "\t\t: jprev@(%g,%g)\n", jprev->vertex.x, jprev->vertex.y);

    assert(jnext!=jprev);
    
    
}