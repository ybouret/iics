#include "rescaler.hpp"


bool Rescaler:: need_to_refine( const Bubble &bubble )
{
    // assume we have a predefined metrics
    const size_t n = bubble.size;
    assert( s.size() == n );
    
    // follow the abscissa
    a_list.empty();
    Real    s_curr   = 0;
    Tracer *p        = bubble.root;
    bool    doRefine = false;
    for( size_t i=n;i>0;--i, p=p->next)
    {
        a_list.append()->s = s_curr;
        const Real ds = p->s;
        if( ds > bubble.lambda )
        {
            // insert supplementary
            doRefine = true;
            a_list.append()->s = s_curr + ds/2;
        }
        s_curr += ds;
    }
    return doRefine;
}


void Rescaler:: refine( Bubble &bubble )
{
    // assume we have a predefined metrics
    const size_t n = bubble.size;
    assert( s.size() == n );
    
    while( need_to_refine(bubble) )
    {
        rebuild(bubble);
    }
    
    
}
