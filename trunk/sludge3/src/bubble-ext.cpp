#include "bubble.hpp"

void Bubble:: auto_contour()
{
    assert(size>=3);
    Tracer *cur = root;
    Tracer *nxt = cur->next;
    std::cerr << "size before =" << size << std::endl;
    const size_t ns = size;
    for(size_t p=ns;p>0;--p)
    {
        const Real dist = cur->dist;
        if( dist > lambda )
        {
            // insert a few points
            const size_t n = size_t( Floor(dist/lambda) ); assert(n>0);
            //std::cerr << "insert " << n << std::endl;
            // with this origin
            const Vertex org = cur->pos;
            Tracer      *ptr = cur;
            const Real   fac = 1.0 / (n+1);
            for(size_t i=1; i <= n; ++i)
            {
                const Vertex v = org + (i*fac) * cur->edge;
                Tracer *tr = new Tracer(v);
                insert_after(ptr,tr);
                assert(owns(ptr));
                assert(owns(tr));
                ptr = tr;
            }
        }
        cur = nxt;
        nxt = nxt->next;
    }
    std::cerr << "size after=" << size << std::endl;
    init_contour();
}
