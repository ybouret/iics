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
    
    
    static inline Real __dist( const Vertex &lhs, const Vertex &rhs )
    {
        return Hypotenuse(rhs.x-lhs.x, rhs.y-lhs.y);
    }
    
    class Adjust
    {
    public:
        const size_t np; //!< #size+1 points (cyclic)
        const Real   L;  //!< perimeter
        vector<Real> t;  //!< 0->L
        matrix<Real> P;  //!< points,derivative,vectors a and b
        
        
        explicit Adjust( const Bubble &b ) :
        np(b.size+1),
        L(0),
        t(np,0),
        P(np,8)
        {
            //------------------------------------------------------------------
            // Fill the coordinates
            //------------------------------------------------------------------
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
            
            {
                ios::ocstream fp("raw.dat",false);
                for(size_t i=1; i <= np; ++i )
                {
                    fp("%g %g %g\n", t[i], P[i][1], P[i][2]);
                }
            }
            
            std::cerr << "L=" << L << std::endl;
            
            //------------------------------------------------------------------
            // Compute the derivatives
            //------------------------------------------------------------------
            for(size_t i=2;i<np;++i)
            {
                const Real   tm   = t[i]   - t[i-1];
                const Real   tp   = t[i+1] - t[i];
                const Vertex Mm   = get_pos(i-1);
                const Vertex Mp   = get_pos(i+1);
                const Vertex M0   = get_pos(i);
                const Vertex dMdt = Derivative(M0, tm, Mm, tp, Mp);
                P[i][3] = dMdt.x;
                P[i][4] = dMdt.y;
            }
            
            {
                const Real   tm = t[np] - t[np-1];
                const Real   tp = t[2]  - t[1];
                const Vertex Mm = get_pos(np-1);
                const Vertex M0 = get_pos(1);
                const Vertex Mp = get_pos(2);
                const Vertex dMdt = Derivative(M0, tm, Mm, tp, Mp);
                P[1][3] = P[np][3] = dMdt.x;
                P[1][4] = P[np][4] = dMdt.y;
            }
            
            {
                ios::ocstream fp("der.dat",false);
                for(size_t i=1; i <= np; ++i )
                {
                    fp("%g %g %g\n", t[i], P[i][1], P[i][2]);
                    fp("%g %g %g\n", t[i], P[i][1]+P[i][3], P[i][2]+P[i][4]);
                    fp("\n");
                }
            }
            
            //------------------------------------------------------------------
            // Compute the order 2 vectors
            //------------------------------------------------------------------
            for(size_t i=1;i<np;++i)
            {
                const Vertex M0 = get_pos(i);
                const Vertex M1 = get_pos(i+1);
                const Real   dt = t[i+1] - t[i];
                const Vertex M0M1(M0,M1);
                
            }
            
        }
        
        inline
        Vertex get_pos(size_t i) const throw()
        {
            return Vertex(P[i][1],P[i][2]);
        }
        
        inline
        Vertex get1( Real u ) const throw()
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
        
        static inline
        Vertex Derivative(const Vertex M0,
                          const Real   tm,
                          const Vertex Mm,
                          const Real   tp,
                          const Vertex Mp
                          ) throw()
        {
            const Vertex M0Mm(M0,Mm);
            const Vertex M0Mp(M0,Mp);
            const Vertex Vm = M0Mm / tm;
            const Vertex Vp = M0Mp / tp;
            const Vertex dV = tm*Vp - tp*Vm;
            return  dV/(tm+tp);
        }
        
#define __CALL(object,ptrToMember)  ((object).*(ptrToMember))
#define __GET(U) ( __CALL(*this,get)(U) )
        
        bool build_ring(Tracer::Ring &ring,
                        const size_t  NP,
                        Vertex (Adjust::*get)(Real) const,
                        const Real lambda
                        )
        {
            ring.auto_delete();
            ring.push_back( new Tracer( __GET(0) ) );
            const Real du   = L/NP;
            
            //------------------------------------------------------------------
            // append points, tracking max distance
            //------------------------------------------------------------------
            for( size_t i=1; i < NP; ++i )
            {
                const Real  u  = i * du;
                Tracer     *tr = new Tracer( __GET(u) );
                ring.push_back(tr);
                const Real  d = __dist(tr->pos,tr->prev->pos);
                if(d>lambda)
                {
                    return false;
                }
                
            }
            
            //------------------------------------------------------------------
            // check last point vs root
            //------------------------------------------------------------------
            const Real  d = __dist(ring.root->pos,ring.root->prev->pos);
            if(d>lambda)
            {
                return false;
            }
            
            return true;
        }
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Adjust);
    };
}







void Bubble:: adjust_contour()
{
    assert(size>=3);
    
    Adjust adjust(*this);
    std::cerr << "Area0=" << area << std::endl;
    
    {
        size_t NP = max_of<Real>(3,adjust.L/lambda);
        Tracer::Ring ring1;
    GENERATE_RING:
        if( ! adjust.build_ring(ring1, NP, & Adjust::get1, lambda) )
        {
            ++NP;
            goto GENERATE_RING;
        }
        swap_with(ring1);
    }
    std::cerr << "Area1=" << __area() << std::endl;
    
    save_dat( "adj.dat" );
    init_contour();
    
}

