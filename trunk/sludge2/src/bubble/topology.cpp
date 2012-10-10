#include "../bubble.hpp"
#include "yocto/code/utils.hpp"

//! simultaneous tridiagonal system
static inline
void vertex_tridiag(const array<Real>   &a,
                    const array<Real>   &b,
                    const array<Real>   &c,
                    const array<Vertex> &r,
                    array<Vertex>       &u,
                    const size_t         n)
{
    const Vertex org(0,0);
    
    assert(n<=a.size());
    assert(n<=b.size());
    assert(n<=c.size());
    assert(n<=r.size());
    assert(n<=u.size());
    assert(Fabs(b[1])>0);
    Vertex         bet;
    vector<Vertex> gam(n,org);
    u[1]=r[1]/(bet.x=bet.y=b[1]);
    for(size_t j=2;j<=n;j++)
    {
        Vertex    &gam_j = gam[j];
        const Real c_jm  = c[j-1];
        gam_j.x = c_jm/bet.x;
        gam_j.y = c_jm/bet.y;
        
        const Real b_j = b[j];
        const Real a_j = a[j];
        
        bet.x = b_j - a_j * gam_j.x;
        bet.y = b_j - a_j * gam_j.y;
        
        if ( Fabs(bet.x) <= 0 || Fabs(bet.y) <= 0 )
            throw exception("Error 2 in tridag");
        
        const Vertex &r_j = r[j];
        Vertex       &u_j = u[j];
        const Vertex &u_jm = u[j-1];
        
        u_j.x = (r_j.x - a_j * u_jm.x) / bet.x;
        u_j.y = (r_j.y - a_j * u_jm.y) / bet.y;
    }
    for(size_t j=(n-1);j>=1;j--)
    {
        Vertex       &u_j    = u[j];
        const Vertex &u_jp   = u[j+1];
        const Vertex &gam_jp = gam[j+1];
        u_j.x -= gam_jp.x * u_jp.x;
        u_j.y -= gam_jp.y * u_jp.y;
    }
    
}


static inline
void vertex_cyclic(const array<Real>   &a,
                   const array<Real>   &b,
                   const array<Real>   &c,
                   const array<Vertex> &r,
                   array<Vertex>       &x,
                   const size_t         n)
{
    assert(n>2);
    assert(n<=a.size());
    assert(n<=b.size());
    assert(n<=c.size());
    assert(n<=r.size());
    assert(n<=x.size());
    const Vertex org(0,0);
    const Real   alpha = c[n];
    const Real   beta  = a[1];
    
    
    vector<Real>   bb(n,0);
    vector<Vertex> u(n,org);
    vector<Vertex> z(n,org);
    
    const Real gamma = -b[1];
    bb[1]=b[1]-gamma;
    bb[n]=b[n]-alpha*beta/gamma;
    for(size_t i=2;i<n;i++)
        bb[i]=b[i];
    vertex_tridiag(a,bb,c,r,x,n);
    u[1].x=u[1].y=gamma;
    u[n].x=u[n].y=alpha;
    for(size_t i=2;i<n;i++)
        u[i].x=u[i].y = 0;
    vertex_tridiag(a,bb,c,u,z,n);
    const Real fact_x= (x[1].x+beta*x[n].x/gamma)/(1.0+z[1].x+beta*z[n].x/gamma);
    const Real fact_y= (x[1].y+beta*x[n].y/gamma)/(1.0+z[1].y+beta*z[n].y/gamma);
    
    for(size_t i=1;i<=n;i++)
    {
        x[i].x -= fact_x*z[i].x;
        x[i].y -= fact_y*z[i].y;
    }
    
}

static inline void __bracket( size_t &klo, size_t &khi, const Real t, const array<Real> &ta)
{
    const size_t n = ta.size();
    assert(t<ta[n]);
    if( (klo<=0) || (klo>=n) || (khi-klo>1) || (t<ta[klo]) || (t>ta[khi]) )
    {
        klo=1;
        khi=n;
        while(khi-klo>1)
        {
            const size_t k=(khi+klo) >> 1;
            if (t<ta[k])
                khi=k;
            else
                klo=k;
        }
    }
    assert(klo>0);
    assert(klo<n);
    assert(1==khi-klo);
    assert(t>=ta[klo]);
    assert(t<=ta[khi]);
}

static inline
Vertex vertex_splint( const Real t, const array<Real> &ta, const array<Vertex> &v, const array<Vertex> &v2 )
{
    static size_t klo = 0;
    static size_t khi = 0;
    assert(ta.size()==v.size());
    assert(ta.size()==v2.size());
    __bracket(klo,khi, t, ta);
    
    const Real h=ta[khi]-ta[klo];
    const Real a=(ta[khi]-t)/h;
    const Real b=(t-ta[klo])/h;
    const Real f= h*h/6.0;
    const Real C=(a*a*a-a)*f;
    const Real D=(b*b*b-b)*f;
    return a*v[klo] + b*v[khi] + C * v2[klo] + D * v2[khi];
}

static inline
Vertex vertex_dsplint1(const Real t, const array<Real> &ta, const array<Vertex> &v, const array<Vertex> &v2 )
{
    static size_t klo = 0;
    static size_t khi = 0;
    assert(ta.size()==v.size());
    assert(ta.size()==v2.size());
    __bracket(klo,khi, t, ta);
    
    const Real h=ta[khi]-ta[klo];
    const Real a=(ta[khi]-t)/h;
    const Real b=(t-ta[klo])/h;
    const Real f=h/6.0;
    const Real Cp=(1.0-3*a*a)*f;
    const Real Dp=(3*b*b-1.0)*f;
    
    return (v[khi] - v[klo])/h + Cp * v2[klo] + Dp * v2[khi];
}

static inline
Vertex vertex_dsplint2(const Real t, const array<Real> &ta, const array<Vertex> &v, const array<Vertex> &v2 )
{
    static size_t klo = 0;
    static size_t khi = 0;
    assert(ta.size()==v.size());
    assert(ta.size()==v2.size());
    __bracket(klo,khi, t, ta);
    
    const Real h=ta[khi]-ta[klo];
    const Real a=(ta[khi]-t)/h;
    const Real b=(t-ta[klo])/h;
    
    return a * v2[klo] + b * v2[khi];
}

void Bubble:: compute_contour()
{
    const Vertex org(0,0);
    
    assert(size>=3);
    const size_t   N = size;
    const size_t   M = N+1;
    vector<Real>   t(M,0);    // curv. parameter
    vector<Vertex> v(M,org);  // absolute verticies
    vector<Vertex> v2(M,org); // for spline coefcicients
    vector<Real>   h(N,0);    // |v[i+1] - v[i]|
    
    //--------------------------------------------------------------------------
    //
    // cumulate both s and contour v
    //
    //--------------------------------------------------------------------------
    Tracer *tracer = root; // starting point
    Real    period = 0;    // current curvilinear absc.
    Vertex  vtx    = org;  // current relative vertex
    
    v[1] = root->vertex;   // PBC ?
    pbc(v[1]);
    for( size_t i=1;i<=N;++i,tracer=tracer->next)
    {
        pbc( tracer->vertex );
        tracer->edge = tracer->next->vertex - tracer->vertex;
        pbc( tracer->edge );
        t[i] = period;
        period += (h[i] = tracer->edge.norm() );
        v[i] = v[1] + vtx;
        vtx += tracer->edge;
    }
    t[M] = period;
    v[M] = v[1];
    
#if 0
    {
        ios::ocstream fp("bv.dat",false);
        for( size_t i=1; i <= M; ++i )
            fp("%g %g %g\n", v[i].x, v[i].y, t[i]);
        
    }
#endif
    
    //--------------------------------------------------------------------------
    //
    // compute spline coefficients
    //
    //--------------------------------------------------------------------------
    vector<Real>   a(N,0);
    vector<Real>   b(N,0);
    vector<Real>   c(N,0);
    vector<Vertex> r(N,org);
    for( size_t i=2; i<N; ++i )
    {
        a[i] =  h[i-1]      / 6.0;
        b[i] = (h[i-1]+h[i])/ 3.0;
        c[i] = h[i]         / 6.0;
        r[i] = (v[i+1] - v[i])/h[i] - (v[i]-v[i-1])/h[i-1];
    }
    
    a[1] = h[N]        / 6.0;
    b[1] = (h[N]+h[1]) / 3.0;
    c[1] = h[1]        / 6.0;
    r[1] = (v[2] - v[1])/h[1] - (v[1]-v[N])/h[N];
    
    a[N] =  h[N-1] / 6.0;
    b[N] = (h[N-1]+h[N]) / 3.0;
    c[N] =  h[1] / 6.0;
    r[N] = (v[M] - v[N])/h[N] - (v[N]-v[N-1])/h[N-1];
    
    
    vertex_cyclic(a, b, c, r, v2, N);
    v2[M] = v2[1];
    
#if 0
    {
        ios::ocstream fp("bs.dat",false);
        for( size_t i=0; i <100; ++i )
        {
            Vertex I = vertex_splint( (period*i)/100, t, v, v2 );
            pbc(I);
            fp("%g %g\n", I.x, I.y);
        }
        fp("%g %g\n", v[1].x, v[1].y);
    }
#endif
    
    //--------------------------------------------------------------------------
    //
    // lazy insertion
    //
    //--------------------------------------------------------------------------
    size_t Ns = max_of<size_t>( size, max_of<size_t>(Ceil( period/lam ),3));
    
TRY_GENERATE:
    {
        //----------------------------
        //-- create temporary tracers
        //----------------------------
        const Real dt = period / Ns;
        Tracers    tmp( cache );
        for( size_t i=0; i < Ns; ++i )
        {
            tmp.append()->vertex = vertex_splint( i * dt, t, v, v2);
        }
        
        //-------------
        //-- check'em
        //-------------
        Tracer *p = tmp.root;
        for(size_t i=tmp.size;i>0;--i,p=p->next)
        {
            p->edge = p->next->vertex - p->vertex;
            p->s    = Sqrt( p->s2 = p->edge.norm2());
            if( p->s > lam)
            {
                ++Ns;
                goto TRY_GENERATE;
            }
        }
        //std::cerr << "generated with Ns=" << Ns << " / " << size << std::endl;
#if 0
        {
            ios::ocstream fp("bg.dat", false);
            Tracer *p = tmp.root;
            for(size_t i=tmp.size;i>0;--i,p=p->next)
            {
                fp("%g %g\n", p->vertex.x, p->vertex.y);
            }
            fp("%g %g\n", p->vertex.x, p->vertex.y);
        }
#endif
        
        // accepted !
        this->swap_with( tmp );
    }
   
    //--------------------------------------------------------------------------
    //
    // pbc and differential properties
    //
    //--------------------------------------------------------------------------
    tracer = root;
    for(size_t i=0; i<size; ++i,tracer=tracer->next)
    {
        tracer->bubble = this;
        pbc(tracer->vertex);
        assert( tracer->edge.norm() <= lam );
        const Real   ti           = (i*period)/size;
        const Vertex tangent      = vertex_dsplint1(ti, t, v, v2);
        const Real   tg_norm2     = tangent.norm2();
        const Real   tg_norm      = Sqrt( tg_norm2 );
        tracer->t                 = (1/tg_norm) * tangent;
        tracer->n.x               = -tracer->t.y;
        tracer->n.y               =  tracer->t.x;
        tracer->angle             =  tracer->n.positive_angle();
        const Vertex acc          =  vertex_dsplint2(ti, t, v, v2);
        tracer->curvature         = (acc * tracer->n) / tg_norm2;
    }
    
    
    //--------------------------------------------------------------------------
    // update area with iso pressure
    //--------------------------------------------------------------------------
    update_area();
    content = area * pressure;
}
