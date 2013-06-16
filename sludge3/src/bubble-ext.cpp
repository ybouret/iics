#include "bubble.hpp"
#include "yocto/exception.hpp"

#if 0
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
#endif

#include "yocto/sequence/vector.hpp"
#include "yocto/math/kernel/matrix.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"
#include <cstring>

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
        
        static const size_t COLS = 8;
        
        explicit Adjust( const Bubble &b ) :
        np(b.size+1),
        L(0),
        t(np,0),
        P(np,COLS)
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
            
#if 0
            {
                ios::ocstream fp("raw.dat",false);
                for(size_t i=1; i <= np; ++i )
                {
                    fp("%g %g %g\n", t[i], P[i][1], P[i][2]);
                }
            }
#endif
            
            //std::cerr << "L=" << L << std::endl;
            
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
            
#if 0
            {
                ios::ocstream fp("der.dat",false);
                for(size_t i=1; i <= np; ++i )
                {
                    fp("%g %g %g\n", t[i], P[i][1], P[i][2]);
                    fp("%g %g %g\n", t[i], P[i][1]+P[i][3], P[i][2]+P[i][4]);
                    fp("\n");
                }
            }
#endif
            //------------------------------------------------------------------
            // Compute the order 2 vectors
            //------------------------------------------------------------------
            for(size_t i=1;i<np;++i)
            {
                const size_t ip = i+1;
                const Vertex M0 = get_pos(i);
                const Vertex M1 = get_pos(ip);
                const Real   dt = t[ip] - t[i];
                const Vertex M0M1(M0,M1);
                const Vertex v0(P[i][3],P[i][4]);      // left  derivative
                const Vertex v1(P[ip][3],P[ip][4]);  // right derivative
                const Vertex U = M0M1 - dt * v0;
                const Vertex V = dt*(v1-v0);
                const Real   dt2 = dt*dt;
                const Real   dt3 = dt2*dt;
                const Vertex a = (3.0*U-V)/dt2;
                const Vertex b = (V-(U+U))/dt3;
                P[i][5] = a.x;
                P[i][6] = a.y;
                P[i][7] = b.x;
                P[i][8] = b.y;
            }
            
            memcpy( &P[np][1], &P[1][1], COLS * sizeof(Real) );
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
            const Vertex M0( P[jlo][1], P[jlo][2]);
            const Vertex M1( P[jhi][1], P[jhi][2]);
            const Real   t0 = t[jlo];
            const Real   t1 = t[jhi];
            const Real   h  = (t1-t0);
            const Vertex M = (t1-u) * M0 + (u-t0) * M1;
            
            return M/h;
        }
        
        inline
        Vertex get2( Real u ) const throw()
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
            const Vertex M0( P[jlo][1], P[jlo][2]);
            const Vertex v0( P[jlo][3], P[jlo][4]);
            const Vertex a(  P[jlo][5], P[jlo][6]);
            const Vertex b(  P[jlo][7], P[jlo][8]);

            const Real   t0 = t[jlo];
            //const Real   t1 = t[jhi];
            const Real   h  = (u-t0);

            return M0 + h * ( v0 + h* (a + h * b) );
        }
        
        inline ~Adjust() throw() {}
        
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
                        const Real lambda,
                        Real      &dmax
                        )
        {
            ring.auto_delete();
            ring.push_back( new Tracer( __GET(0) ) );
            const Real du   = L/NP;
            
            //------------------------------------------------------------------
            // append points, tracking max distance
            //------------------------------------------------------------------
            dmax = 0;
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
                if(d>dmax) dmax=d;
                
            }
            
            //------------------------------------------------------------------
            // check last point vs root
            //------------------------------------------------------------------
            const Real  d = __dist(ring.root->pos,ring.root->prev->pos);
            if(d>lambda)
            {
                return false;
            }
            if(d>dmax)
                dmax=d;
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
    const Real A0 = area;    
    {
        size_t NP = max_of<Real>(3,adjust.L/lambda);
        Tracer::Ring ring;
        Real         dmax = 0;
        
        //======================================================================
        //
        // Subdivide
        //
        //======================================================================
    GENERATE_RING:
        if( !adjust.build_ring(ring, NP, & Adjust::get2, lambda, dmax))
        {
            ++NP;
            goto GENERATE_RING;
        }
        
        //======================================================================
        //
        // New Area => expansion factopr
        //
        //======================================================================

        const Real A1    = ring.__area();       
        const Real ratio = Sqrt(A0/A1);
        const Real expanded = ratio * dmax;
        
        //======================================================================
        //
        // Finalize
        //
        //======================================================================

        if( expanded >= lambda )
        {
            ++NP;
            goto GENERATE_RING;
        }
        
        swap_with(ring);
        ring.auto_delete();
        assert(size==NP);
        
        //======================================================================
        //
        // Expand
        //
        //======================================================================
        Tracer *tr = root;
        for(size_t i=size;i>0;--i,tr=tr->next)
        {
            const Vertex dr(G,tr->pos);
            tr->pos = G + ratio * dr;
        }
    }
    
    init_contour();
    
}

