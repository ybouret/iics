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
idelta(0),
idelta2(0)
{
}


Vertex Arc:: operator()( const Real mu ) const throw()
{
    return r0 + mu * a + (mu*mu) * b + (mu*mu*mu) * c;
}

#include "yocto/code/ipower.hpp"



#define ARC_NVAR 4

void Arc:: load( const array<Real> &U ) const
{
    assert( U.size() == ARC_NVAR );
    a.x = U[1];
    a.y = U[2];
    b.x = U[3];
    b.y = U[4];
    c.x = delta_r.x - (a.x+b.x);
    c.y = delta_r.y - (a.y+b.y);
}

void Arc:: init(array<Real> &U) const
{
    delta_r = r1 - r0;
    delta   = delta_r.norm();
    idelta  = 1/delta;
    idelta2 = 1/(delta*delta);
    
    const Real scale = 0.5 * delta;
    a = scale * t0;
    b = 3.0*delta_r - scale * ( t1+2.0*t0);
    
    U[1] = a.x;
    U[2] = a.y;
    U[3] = b.x;
    U[4] = b.y;
    load(U);
    std::cerr << "init=" << U << std::endl;
}


void Arc:: func( array<Real> &F, const array<Real> &U) const
{
    load(U);
    
    const Vertex dr0   = a;
    const Real   dr0n  = dr0.norm();
    const Real   dr0n3 = dr0n * dr0n * dr0n;
    
    const Vertex dr1   = a+2.0*b+3.0*c;
    const Real   dr1n  = dr1.norm();
    const Real   dr1n3 = dr1n * dr1n * dr1n;
    
    F[1] = 1 - (dr0*t0)/dr0n;
    F[2] = 1 - (dr1*t1)/dr1n;
    F[3] = delta * ( 0.5*C0 - Vertex::det(a,b)/dr0n3 );
    F[4] = delta * ( 0.5*C1 - (Vertex::det(a,b)+3*Vertex::det(b,c)+3*Vertex::det(a,c))/dr1n3);
}

void Arc:: fjac( matrix<Real> &J, const array<Real> &U ) const
{
    load(U);
    J.ldz();
    
    const Vertex dr0   = a;
    const Real   dr0n  = dr0.norm();
    const Real   dr0n2 = dr0n * dr0n;
    const Real   dr0n3 = dr0n * dr0n2;
    const Real   as0   = dr0 * t0;
    
    const Vertex dr1   = a+2.0*b+3.0*c;
    const Real   dr1n  = dr1.norm();
    const Real   dr1n2 = dr1n * dr1n;
    const Real   dr1n3 = dr1n * dr1n2;
    const Real   as1   = dr1 * t1;
    std::cerr << "[dr0n=" << dr0n << ", dr1n=" << dr1n << " ]" << std::endl;
    
    {
        array<Real> &Q = J[1];
        Q[1] = ((a.x*as0)/dr0n2 - t0.x)/dr0n;
        Q[2] = ((a.y*as0)/dr0n2 - t0.y)/dr0n;
    }
    
    {
        array<Real> &Q = J[2];
        const Real  gx = (t1.x - (dr1.x * as1) /dr1n2)/dr1n;
        const Real  gy = (t1.y - (dr1.y * as1) /dr1n2)/dr1n;
        Q[1] = 2*gx;
        Q[2] = 2*gy;
        Q[3] = gx;
        Q[4] = gy;
    }
    
    {
        array<Real> &Q = J[3];
        const Real fac   = 3* Vertex::det(a,b)/dr0n2;
        const Real scale = delta / dr0n3;
        Q[1] = scale * ( dr0.x * fac - b.y );
        Q[2] = scale * ( dr0.y * fac + b.x );
        Q[3] = scale * (  a.y );
        Q[4] = scale * ( -a.x );
    }
    
    {
        array<Real> &Q   = J[4];
        const Real fac   = 3 * (  Vertex::det(a,b) + 3*Vertex::det(a,delta_r) + 3*Vertex::det(b,delta_r) )/dr1n2;
        const Real scale = delta / dr1n3;
        const Real gx    = dr1.x * fac;
        const Real gy    = dr1.y * fac;
        
        Q[1] = - scale * (  3*delta_r.y+b.y + 2*gx);
        Q[2] = - scale * ( -3*delta_r.x-b.x + 2*gy);
        Q[3] = - scale * (  3*delta_r.y-a.y +   gx);
        Q[4] = - scale * ( -3*delta_r.x+a.x +   gy);
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

