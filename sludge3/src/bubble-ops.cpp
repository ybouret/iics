#include "bubble.hpp"

void Bubble:: init_contour() throw()
{
    assert(size>=3);
    Tracer *tr = root;
    for( size_t i=size;i>0;--i,tr=tr->next)
    {
        tr->edge = tr->next->pos - tr->pos;
        tr->dist = tr->edge.norm();
        std::cerr << "\td=" << tr->dist << std::endl;
    }
}

#include "yocto/math/dat/spline.hpp"
#include "yocto/ios/ocstream.hpp"

namespace {
    
    class Spline
    {
    public:
        const size_t ns;
        vector<Real> t;
        matrix<Real> P;
        matrix<Real> Q;
        const Real   width;
        
        inline Spline( const Tracer::Ring &ring ) :
        ns( ring.size + 1),
        t(ns,0),
        P(2,ns),
        Q(2,ns),
        width(0)
        {
            array<Real>  &X  = P[1];
            array<Real>  &Y  = P[2];
            const Tracer *tr = ring.root;
            
            t[1] = 0;
            X[1] = tr->pos.x;
            Y[1] = tr->pos.y;
            tr = tr->next;
            
            for(size_t i=ring.size,j=2;i>0;--i,++j,tr=tr->next)
            {
                t[j] = t[j-1] + tr->dist;
                X[j] = tr->pos.x;
                Y[j] = tr->pos.y;
            }
            
            (Real &)width = t[ns];
            spline2D<Real>::compute(spline_periodic, t, P, Q);
        }
        
        inline ~Spline() throw()
        {
        }
        
        inline Vertex operator()( Real u ) const throw()
        {
            return spline2D<Real>::eval(u, t, P, Q, &width);
        }
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Spline);
    };
}

#include "yocto/code/utils.hpp"

static inline
bool are_valid( const Tracer *p, const Tracer *q, const Real lam)
{
    return Hypotenuse(p->pos.x-q->pos.x,
                      p->pos.y-q->pos.y) <= lam;
}

void Bubble:: auto_contour()
{
    assert(size>=3);
    const size_t ns = size+1;
    
    //==========================================================================
    //
    // compute the spline2D
    //
    //==========================================================================
    Spline S(*this);
    
    {
        ios::ocstream fp( "poly.dat", false );
        for(size_t i=1; i <= ns; ++i )
        {
            fp("%g %g %g\n", S.t[i], S.P[1][i], S.P[2][i]);
        }
    }
    
    {
        ios::ocstream fp( "spoly.dat", false );
        const size_t NS = 200;
        for( size_t i=0; i < NS; ++i )
        {
            const Real   u = (i*S.width) / NS;
            const Vertex v = S(u);
            fp("%g %g %g\n", u, v.x, v.y);
        }
    }
    
    
    //==========================================================================
    //
    // Remap by lazy insertion
    //
    //==========================================================================
    
    //--------------------------------------------------------------------------
    // initialize the ring
    //--------------------------------------------------------------------------
    const Vertex org  = root->pos;
    Tracer::Ring ring;
    
    
    //--------------------------------------------------------------------------
    // fill with distance control
    //--------------------------------------------------------------------------
    size_t N = max_of<size_t>(3,S.width / lambda);
GENERATE:
    {
        std::cerr << "\t\tmapping with " << N << " points" << std::endl;
        ring.push_back( new Tracer(org) );
        const Real du = S.width / N;
        for(size_t i=1;i<N;++i)
        {
            Tracer *tr = new Tracer( S(du*i) );
            ring.push_back(tr);
            assert(0!=tr->prev);
            if( ! are_valid(tr, tr->prev, lambda) )
            {
                ++N;
                ring.auto_delete();
                goto GENERATE;
            }
        }
        assert(ring.size==N);
        if( ! are_valid(ring.root, ring.root->prev, lambda) )
        {
            ++N;
            ring.auto_delete();
            goto GENERATE;
        }
    }
    
    //--------------------------------------------------------------------------
    // winner
    //--------------------------------------------------------------------------
    ring.swap_with(*this);
    
    
    
}

