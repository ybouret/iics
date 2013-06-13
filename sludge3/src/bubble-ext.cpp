#include "bubble.hpp"
#include "yocto/exception.hpp"

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



bool Bubble:: reduce1(const Real mu)
{
    if(size>=4)
    {
        Tracer *t1 = root;
        Tracer *t2 = t1->next;
        Tracer *t3 = t2->next;
        Tracer *t4 = t3->next;
        
        for(size_t i=size;i>0;--i)
        {
            const Vertex &M2 = t2->pos;
            const Vertex &M3 = t3->pos;
            const Vertex  M2M3(M2,M3);
            const Real    d23 = M2M3.norm();
            std::cerr << "d23=" << d23 << "/mu=" << mu << std::endl;
            if(d23<=mu)
            {
                const Vertex t = M2M3/d23;
                const Vertex v(-t.y,t.x);
                const Vertex &M1 = t1->pos;
                const Vertex &M4 = t4->pos;
                const Vertex  M1M4(M1,M4);
                const Real   A0 = Vertex::det(M1,M2) + Vertex::det(M2,M3) + Vertex::det(M3,M4);
                const Vertex I  = 0.5 * ( M2+M3 );
                const Real   den = Vertex::det(v,M1M4);
                if(Fabs(den)<=numeric<Real>::minimum)
                    throw exception("Ill-Formed Polygon for Reduction");
                
                const Real   alpha = (A0 + Vertex::det(M1M4,I))/den;
                const Vertex Q     = I + alpha * v;
                std::cerr << "alpha=" << alpha << ", v=" << v << ": " << I << " -> " << Q << std::endl;

                std::cerr << "A0=" << A0 << std::endl;
                std::cerr << "A1=" << Vertex::det(M1,Q) + Vertex::det(Q,M4) << std::endl;
                Tracer *tr = new Tracer(Q);
                std::cerr << "Area0=" << __area() << std::endl;
                delete unlink(t2);
                delete unlink(t3);
                insert_after(t1, tr);
                std::cerr << "Area1=" << __area() << std::endl;
                return true;
            }
            t1=t2;
            t2=t3;
            t3=t4;
            t4=t4->next;
        }
    }
    return false;
}

void Bubble:: reduce( const Real mu )
{
    while( reduce1(mu) )
        ;
    init_contour();
}
