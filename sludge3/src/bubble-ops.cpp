#include "bubble.hpp"

Real Bubble:: __area() const throw()
{
    Real ans = 0;
    Tracer *tr = root;
    for( size_t i=size;i>0;--i,tr=tr->next)
    {
        const Vertex &p = tr->pos;
        const Vertex &q = tr->next->pos;
        ans += p.x * q.y - p.y * q.x;
    }
    return Fabs(ans)/2;
}


void Bubble:: init_contour() throw()
{
    assert(size>=3);
    G.ldz();
    area = 0;
    Tracer *tr = root;
    for( size_t i=size;i>0;--i,tr=tr->next)
    {
        const Vertex &p = tr->pos;
        G += p;
        const Vertex &q = tr->next->pos;
        tr->edge = q - p;
        tr->dist = tr->edge.norm();
        area += p.x * q.y - p.y * q.x;
    }
    G.x /= size;
    G.y /= size;
    area = Fabs(area)/2;
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
            tr   = tr->next;
            
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
Real __distance( const Tracer *p, const Tracer *q ) throw()
{
    assert(p!=q);
    const Real dx = p->pos.x - q->pos.x;
    const Real dy = p->pos.y - q->pos.y;
    return Sqrt(dx*dx+dy*dy);
}

static inline
bool __are_valid( const Tracer *p, const Tracer *q, const Real lam)
{
    return __distance(p,q) <= lam;
}

void Bubble:: auto_contour_spline()
{
    assert(size>=3);
    assert(area>0);
    
    //==========================================================================
    //
    // compute the spline2D
    //
    //==========================================================================
    Spline S(*this);
    
#if 1
    const size_t ns = size+1;
    {
        ios::ocstream fp( "poly.dat", false );
        for(size_t i=1; i <= ns; ++i )
        {
            fp("%g %g %g\n", S.t[i], S.P[1][i], S.P[2][i]);
        }
    }
    
    {
        ios::ocstream fp( "spoly.dat", false );
        const size_t NS = 2048;
        for( size_t i=0; i < NS; ++i )
        {
            const Real   u = (i*S.width) / NS;
            const Vertex v = S(u);
            fp("%g %g %g\n", u, v.x, v.y);
        }
    }
#endif
    
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
    std::cerr << "Spline width=" << S.width << std::endl;
    size_t N = max_of<size_t>(3,S.width / lambda);
GENERATE:
    {
        std::cerr << "\t\tmapping " << lambda << " with " << N << " points" << std::endl;
        
        //----------------------------------------------------------------------
        // sampling rate
        //----------------------------------------------------------------------
        const Real du = S.width / N;
        
        //----------------------------------------------------------------------
        // first tracer at origin
        //----------------------------------------------------------------------
        ring.push_back( new Tracer(org) );
        
        //----------------------------------------------------------------------
        // other tracers on spline
        //----------------------------------------------------------------------
        const Tracer *prv = ring.root;
        for(size_t i=1;i<N;++i)
        {
            const Real   u = du * i;
            const Vertex v = S(u);
            Tracer *tr = new Tracer( v );
            ring.push_back(tr);
            assert(0!=tr->prev);
            assert(prv==tr->prev);
            if( ! __are_valid(tr, tr->prev, lambda) )
            {
                std::cerr << "core: invalid " << tr->pos << "/" << tr->prev->pos << ": dist=" << __distance(tr,tr->prev) << std::endl;
                ++N;
                ring.auto_delete();
                goto GENERATE;
            }
            prv = tr;
        }
        assert(ring.size==N);
        assert(ring.root->prev == prv);
        if( ! __are_valid(ring.root, ring.root->prev, lambda) )
        {
            std::cerr << "core: invalid " << ring.root->pos << "/" << ring.root->prev->pos << std::endl;

            ++N;
            ring.auto_delete();
            goto GENERATE;
        }
    }
    
    //--------------------------------------------------------------------------
    // winner
    //--------------------------------------------------------------------------
    ring.swap_with(*this);
    init_contour();
    
}


void Bubble:: compute_curvatures()
{
    assert(size>=3);
    Tracer *tr = root;
    for(size_t i=size;i>0;--i,tr=tr->next)
    {
        tr->compute_order1();
    }
    
    assert(root==tr);
    for(size_t i=size;i>0;--i,tr=tr->next)
    {
        tr->compute_order2();
    }

}


void Bubble:: regularize()
{
    auto_contour();
    compute_curvatures();
}
