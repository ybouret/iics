#include "arc.hpp"


Arc:: ~Arc() throw() {}

Arc:: Arc() throw() :
Func( this, & Arc::func),
Grad( this, & Arc::grad),
a(),
b(),
c(),
r0(),
t0(),
C0(0),
r1(),
t1(),
C1(0),
delta_r(),
delta(0),
delta2(0)
{
}


Vertex Arc:: operator()( const Real mu ) const throw()
{
    return r0 + mu * a + (mu*mu) * b + (mu*mu*mu) * c;
}

#include "yocto/code/ipower.hpp"

Real Arc:: func( const array<Real> &U ) const
{
    Real ans = 0;
    load(U);
    
    {
        const Vertex dpos = ( (a+b+c) - delta_r)/delta;
        ans += 0.5 * dpos.norm2();
    }
    const Vertex dr0  = a;
    const Real   dr0n = dr0.norm();
    {
        ans += (1.0 - (t0*dr0)/dr0n);
    }
    
    const Vertex dr1  = a+2.0*b+3.0*c;
    const Real   dr1n = dr1.norm();
    {
        
        ans += (1.0 - (t1*dr1)/dr1n);
    } 
    
    {
        const Real arg = Vertex::det(a,b)/ipower(dr0n,3) - 0.5*C0;
        ans += 0.5 * delta2 * arg * arg;
    }
    
    {
        const Real arg = (Vertex::det(a,b)+3*Vertex::det(a,c)+3*Vertex::det(b,c))/ipower(dr1n,3) - 0.5*C1;
        ans += 0.5 * delta2 * arg * arg;
    }
    return ans;
}


void Arc:: grad( array<Real> &G, const array<Real> &U) const
{
    assert( G.size() == U.size() );
    for( size_t i=G.size();i>0;--i) G[i] = 0;
    
    
    {
        const Vertex g = ( (a+b+c) - delta_r)/delta2;
        G[1] += g.x;
        G[2] += g.y;
        G[3] += g.x;
        G[4] += g.y;
        G[5] += g.x;
        G[6] += g.y;
    }
    
    const Vertex dr0   = a;
    const Real   dr0n  = dr0.norm();
    const Real   dr0n3 = dr0n*dr0n*dr0n;
    {
        const Real   gs    = dr0 * t0;
        const Real   gx    = dr0.x * gs / dr0n3 - t0.x/dr0n;
        const Real   gy    = dr0.y * gs / dr0n3 - t0.y/dr0n;
        G[1] += gx;
        G[2] += gy;
    }
    
    const Vertex dr1   = a+2.0*b+3.0*c;
    const Real   dr1n  = dr1.norm();
    const Real   dr1n3 = dr1n*dr1n*dr1n;
    
    {
        const Real   gs    = dr1 * t1;
        const Real   gx    = dr1.x * gs / dr1n3 - t1.x/dr1n;
        const Real   gy    = dr1.y * gs / dr1n3 - t1.y/dr1n;
        G[1] += gx;
        G[2] += gy;
        G[3] += 2*gx;
        G[4] += 2*gy;
        G[5] += 3*gx;
        G[6] += 3*gy;
    }
    
    {
        const Real detAB = Vertex::det(a,b);
        const Real arg   = delta2 * (detAB/dr0n3 - 0.5*C0);
        const Real dr0n5 = dr0n3 * dr0n * dr0n;
        G[1] += arg * ( b.y / dr0n3 - 3*a.x*detAB/dr0n5);
        G[2] += arg * (-b.x / dr0n3 - 3*a.y*detAB/dr0n5);
        G[3] += arg * (-a.y / dr0n3 );
        G[4] += arg * ( a.x / dr0n3  );
    }

    {
        const Real detAB = Vertex::det(a,b);
        const Real detBC = Vertex::det(b,c);
        const Real detAC = Vertex::det(a,c);
        const Real sdet  = detAB + 3*(detBC+detAC);
        const Real dr1n5 = dr1n3 * dr1n * dr1n;
        const Real arg   = delta2 * ( sdet/dr1n3 - 0.5*C1);
        const Real fac_x = 3*dr1.x*sdet/dr1n5;
        const Real fac_y = 3*dr1.y*sdet/dr1n5;

        G[1] += arg * ( ( 3*c.y +   b.y )/dr1n3 -     fac_x);
        G[2] += arg * ( (-3*c.x -   b.x )/dr1n3 -     fac_y);
        
        G[3] += arg * ( ( 3*c.y -   a.y )/dr1n3 - 2 * fac_x);
        G[4] += arg * ( (   a.x - 3*c.x )/dr1n3 - 2 * fac_y);

        G[5] += arg * ( (-3*b.y - 3*a.y )/dr1n3 - 3 * fac_x);
        G[6] += arg * ( ( 3*b.x + 3*a.x )/dr1n3 - 3 * fac_y);
        
    }
    
}


#define ARC_NVAR 6

void Arc:: load( const array<Real> &U ) const
{
    assert( U.size() == ARC_NVAR );
    a.x = U[1];
    a.y = U[2];
    b.x = U[3];
    b.y = U[4];
    c.x = U[5];
    c.y = U[6];
}

void Arc:: init(array<Real> &U) const
{
    delta_r = r1 - r0;
    delta2  = delta_r.norm2();
    delta   = Sqrt(delta2);
    
    const Real scale = 0.1 * delta;
    a = scale * t0;
    
#if 0
    // 3 params
    const Vertex v1 = delta_r - a;
    const Vertex v2 = scale * t1 - a;
    
    c = v2 - 2.0*v1;
    b = 3.0*v1 - v2;
#else
    // 2 params
    c.ldz();
    b = 0.5 * (scale*t1 - a );
#endif
    
    U[1] = a.x;
    U[2] = a.y;
    U[3] = b.x;
    U[4] = b.y;
    U[5] = c.x;
    U[6] = c.y;
    std::cerr << "init=" << U << std::endl;
}


bool Arc:: disp( const array<Real> &U ) const
{
    std::cerr << "F(" << U << ")=" << Func(U) << std::endl;
    
    vector<Real> G(U.size(),0);
    Grad(G,U);
    std::cerr << "Grad=" << G << std::endl;
    
    return true;
}


ArcSolver:: ~ArcSolver() throw() {}

ArcSolver:: ArcSolver() :
U(ARC_NVAR,0.0)
{
    
}

#include "yocto/math/kernel/algebra.hpp"
#include "yocto/math/opt/cgrad.hpp"

void ArcSolver:: compute(const Arc &arc)
{
    //--------------------------------------------------------------------------
    // initialize
    //--------------------------------------------------------------------------
    arc.init(U);
    
    cgrad<Real>::callback cb( &arc, &Arc::disp);
    cgrad<Real>::optimize(arc.Func, arc.Grad, U, 1e-5, &cb);
    
    arc.load(U);
    std::cerr << "a=" << arc.a << std::endl;
    std::cerr << "b=" << arc.b << std::endl;
    std::cerr << "c=" << arc.c << std::endl;
}

