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



#define ARC_NVAR 4

void Arc:: load( const array<Real> &U ) const
{
    assert( U.size() == ARC_NVAR );
    a.x = U[1];
    a.y = U[2];
    b.x = U[3];
    b.y = U[4];
    
    //c.x = U[5];
    //c.y = U[6];
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
    //U[5] = c.x;
    //U[6] = c.y;
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
    }
    
    {
        array<Real> &Q = J[2];
        Q[1] = 0;
        Q[2] = idelta;
        Q[3] = 0;
        Q[4] = idelta;
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
        //Q[5] = 3 * Q[1];
        //Q[6] = 3 * Q[2];
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

