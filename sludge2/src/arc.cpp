#include "arc.hpp"


Arc:: ~Arc() throw() {}

Arc:: Arc() throw() :
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
    
    const Real scale = 0.5 * delta;
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
    
    a=delta_r;
    b.ldz();
    c.ldz();
    
    U[1] = a.x;
    U[2] = a.y;
    U[3] = b.x;
    U[4] = b.y;
    U[5] = c.x;
    U[6] = c.y;
    std::cerr << "init=" << U << std::endl;
}


void Arc:: func( array<Real> &F, const array<Real> &U) const
{
    load(U);
    F[1] = (a.x+b.x+c.x - delta_r.x)/delta;
    F[2] = (a.y+b.y+c.y - delta_r.y)/delta;
    
    const Vertex dr0   = a;
    const Real   dr0n  = dr0.norm();
    F[3] = 1 - (dr0*t0)/dr0n;
    
    
    const Vertex dr1 = a+2.0*b+3.0*c;
    const Real   dr1n = dr1.norm();
    F[4] = 1 - (dr1*t1)/dr1n;
    
    //std::cerr << "dr0=" << dr0 << std::endl;
    //std::cerr << "dr1=" << dr1 << std::endl;
    
    F[5] = delta * ( C0/2 - Vertex::det(a,b) / (dr0n*dr0n*dr0n) );
    
    F[6] = delta * ( C1/2 - (Vertex::det(a,b)+3*Vertex::det(b,c)+3*Vertex::det(a,c))/(dr1n*dr1n*dr1n));
}

void Arc:: fjac( matrix<Real> &J, const array<Real> &U ) const
{
    load(U);
    J.ldz();
    const Real idelta = 1/delta;
    
    {
        array<Real> &Q = J[1];
        Q[1] = idelta;
        Q[2] = 0;
        Q[3] = idelta;
        Q[4] = 0;
        Q[5] = idelta;
        Q[6] = 0;
    }
    
    {
        array<Real> &Q = J[2];
        Q[1] = 0;
        Q[2] = idelta;
        Q[3] = 0;
        Q[4] = idelta;
        Q[5] = 0;
        Q[6] = idelta;
    }
    
    const Vertex dr0   = a;
    const Real   dr0n  = dr0.norm();
    const Real   dr0n3 = dr0n * dr0n * dr0n;
    const Real   as0   = dr0 * t0;
    {
        array<Real> &Q = J[3];
        Q[1] = dr0.x * as0 / dr0n3 - t0.x / dr0n;
        Q[2] = dr0.y * as0 / dr0n3 - t0.y / dr0n;
    }
    
    const Vertex dr1   = a+2.0*b+3.0*c;
    const Real   dr1n  = dr1.norm();
    const Real   dr1n3 = dr1n * dr1n * dr1n;
    const Real   as1   = dr1 * t1;
    {
        array<Real> &Q = J[4];
        Q[1] = dr1.x * as1 / dr1n3 - t1.x / dr1n;
        Q[2] = dr1.y * as1 / dr1n3 - t1.y / dr1n;
        Q[3] = 2 * Q[1];
        Q[4] = 2 * Q[2];
        Q[5] = 3 * Q[1];
        Q[6] = 3 * Q[2];
    }
    
    {
        array<Real> &Q = J[5];
        const Real f5 = 3*Vertex::det(a,b)/(dr0n3*dr0n*dr0n);
        Q[1] =  delta*( a.x * f5 - b.y / dr0n3);
        Q[2] =  delta*( a.y * f5 + b.x / dr0n3);
        Q[3] =  delta*a.y/dr0n3;
        Q[4] = -delta*a.x/dr0n3;
    }
    
    std::cerr << "[dr0n=" << dr0n << ", dr1n=" << dr1n << " ]" << std::endl;
    {
        array<Real> &Q = J[6];
        const Real f6 = 3*(Vertex::det(a,b)+3*Vertex::det(b,c)+3*Vertex::det(a,c))/(dr1n3*dr1n*dr1n);
        const Real gx = dr1.x * f6;
        const Real gy = dr1.y * f6;
        Q[1] = delta * (   gx - ( 3*c.y +   b.y)/dr1n3);
        Q[2] = delta * (   gy - (-3*c.x -   b.x)/dr1n3);
        Q[3] = delta * ( 2*gx - ( 3*c.y -   a.y)/dr1n3);
        Q[4] = delta * ( 2*gy - (   a.x - 3*c.x)/dr1n3);
        Q[5] = delta * ( 3*gx - (-3*b.y - 3*a.y)/dr1n3);
        Q[6] = delta * ( 3*gy - ( 3*b.x + 3*a.x)/dr1n3);
    }
}


ArcSolver:: ~ArcSolver() throw() {}

ArcSolver:: ArcSolver() :
U(ARC_NVAR,0.0)
{
    
}

#include "yocto/math/kernel/algebra.hpp"
#include "yocto/math/fcn/newton.hpp"

void ArcSolver:: compute(const Arc &arc)
{
    //--------------------------------------------------------------------------
    // initialize
    //--------------------------------------------------------------------------
    arc.init(U);
    
    Newton<Real>::Function Fn( &arc, & Arc::func );
    Newton<Real>::Jacobian Jn( &arc, & Arc::fjac );
    
    Newton<Real>::solve(Fn, Jn, U, 1e-7);
    
    arc.load(U);
    std::cerr << "a=" << arc.a << std::endl;
    std::cerr << "b=" << arc.b << std::endl;
    std::cerr << "c=" << arc.c << std::endl;
}

