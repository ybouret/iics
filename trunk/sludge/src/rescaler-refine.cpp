#include "rescaler.hpp"

#include "yocto/ios/ocstream.hpp"

bool Rescaler:: need_to_refine( const Bubble &bubble )
{
    //--------------------------------------------------------------------------
    // assume we have a predefined metrics
    //--------------------------------------------------------------------------

    const size_t n = bubble.size;
    assert( s.size()  == n );
    assert( ax.size() == n );
    assert( ay.size() == n );
    assert( theta.size() == n );
    assert(period>0);
    assert(bubble.area>0);
    
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
    
    if( doRefine )
    {
        std::cerr << "s_old=" << s << "/" << period << std::endl;
        std::cerr << "s_new=[";
        for( const abscissa *a=a_list.head;a;a=a->next)
        {
            std::cerr << " " << a->s;
        }
        std::cerr << "]" << std::endl;
    }
    
    return doRefine;
}


void Rescaler:: refine( Bubble &bubble )
{
    static int fid = 0;

    //--------------------------------------------------------------------------
    // assume we have a predefined metrics
    //--------------------------------------------------------------------------
    const size_t n = bubble.size;
    assert( s.size() == n );
    
    while( need_to_refine(bubble) )
    {
        rebuild(bubble);
        ios::ocstream fp( vformat("refine%d.dat",fid++), false);
        const Tracer *p = bubble.root;
        for(size_t i=bubble.size;i>0;--i,p=p->next )
        {
            fp("%g %g\n", p->vertex.x, p->vertex.y);
        }
        fp("%g %g\n", p->vertex.x, p->vertex.y);

    }
    
    
}
