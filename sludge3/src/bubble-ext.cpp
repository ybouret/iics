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


#include "yocto/sequence/vector.hpp"
#include "yocto/math/kernel/matrix.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"

namespace {
    
    class Adjust
    {
    public:
        const size_t np; //!< #size+1 points (cyclic)
        const Real   L;  //!< perimeter
        vector<Real> t;  //!< 0->L
        matrix<Real> P;  //!< points
        
        
        explicit Adjust( const Bubble &b ) :
        np(b.size+1),
        L(0),
        t(np,0),
        P(np,4)
        {
            Tracer *curr = b.root;
            t[1] = 0;
            P[1][1] = curr->pos.x;
            P[1][2] = curr->pos.y;
            for( size_t i=2; i <= np; ++i )
            {
                assert(curr->dist>0);
                t[i] = t[i-1] + curr->dist; //!< distance to next
                curr = curr->next;
                P[i][1] = curr->pos.x;
                P[i][2] = curr->pos.y;
            }

            (Real&)L = t[np];
            
            ios::ocstream fp("raw.dat",false);
            for(size_t i=1; i <= np; ++i )
            {
                fp("%g %g %g\n", t[i], P[i][1], P[i][2]);
            }

            std::cerr << "L=" << L << std::endl;
        }
        
        Vertex operator()( Real u )
        {
            if(u<=0||u>=L)
                return Vertex(P[1][1],P[1][2]);
            
            size_t jlo = 1;
            size_t jhi = np;
            while(jhi-jlo>1)
            {
                const size_t mid = (jlo+jhi)>>1;
                const Real   tmp = t[mid];
                if( tmp > u )
                {
                    jhi = mid;
                }
                else
                {
                    jlo = mid;
                }
            }
            const Vertex v0( P[jlo][1], P[jlo][2]);
            const Vertex v1( P[jhi][1], P[jhi][2]);
            const Real   t0 = t[jlo];
            const Real   t1 = t[jhi];
            const Real   h  = (t1-t0);
            const Vertex v = (t1-u) * v0 + (u-t0) * v1;
            
            return v/h;
        }
        
        ~Adjust() throw() {}
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Adjust);
    };
}

static inline Real __dist( const Vertex &lhs, const Vertex &rhs )
{
    return Hypotenuse(rhs.x-lhs.x, rhs.y-lhs.y);
}

void Bubble:: adjust_contour()
{
    assert(size>=3);
   
    Adjust adjust(*this);
    std::cerr << "Area0=" << area << std::endl;
    
    size_t NP = max_of<Real>(3,adjust.L/lambda);
    {
        Tracer::Ring ring;
        ring.push_back( new Tracer( adjust(0) ) );
        const Real du = adjust.L/NP;
        Real       dmax = 0;
        for( size_t i=1; i < NP; ++i )
        {
            const Real u  = i * du;
            Tracer     *tr = new Tracer( adjust(u) );
            ring.push_back(tr);
            dmax  = max_of(dmax,__dist(tr->pos,tr->prev->pos));
        }
        swap_with(ring);
    }
    std::cerr << "Area1=" << __area() << std::endl;

    save_dat( "adj.dat" );
    init_contour();
    
}

