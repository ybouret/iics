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
delta_r()
{
}


Vertex Arc:: operator()( const Real mu ) const throw()
{
    return r0 + mu * a + (mu*mu) * b + (mu*mu*mu) * c;
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
    const Real scale = 0.1 * delta_r.norm();
    a = scale * t0;
    const Vertex v1 = delta_r - a;
    const Vertex v2 = scale * t1 - a;
    
    c = v2 - 2.0*v1;
    b = v1 - c;
    
    U[1] = a.x;
    U[2] = a.y;
    U[3] = b.x;
    U[4] = b.y;
    U[5] = c.x;
    U[6] = c.y;
    std::cerr << "init=" << U << std::endl;
}

#include "yocto/code/ipower.hpp"

void Arc:: estimate( matrix<Real> &P, array<Real> &F, const array<Real> &U ) const
{
    load(U);
    
    {
        F[1] = a.x+b.x+c.x - delta_r.x;
        array<Real> &J = P[1];
        J[1] = 1;
        J[2] = 0;
        J[3] = 1;
        J[4] = 0;
        J[5] = 1;
        J[6] = 0;
    }
    
    {
        F[2] = a.y+b.y+c.y - delta_r.y;
        array<Real> &J = P[2];
        J[1] = 0;
        J[2] = 1;
        J[3] = 0;
        J[4] = 1;
        J[5] = 0;
        J[6] = 1;
    }
    
    const Vertex dr0      = a;
    const Real   dr0_norm = dr0.norm();
    {
        F[3] = dr0 * t0 - dr0_norm;
        array<Real> &J = P[3];
        J[1] = t0.x - dr0.x / dr0_norm;
        J[2] = t0.y - dr0.y / dr0_norm;
        J[3] = J[4] = J[5] = J[6] = 0;
    }
    
    const Vertex dr1 = a + (b+b) + (c+c+c);
    const Real   dr1_norm = dr1.norm();
    {
        F[4] = dr1 * t1 - dr1_norm;
        array<Real> &J = P[4];
        J[1] = t1.x - dr1.x / dr1_norm;
        J[2] = t1.y - dr1.y / dr1_norm;
        J[3] = 2 * J[1];
        J[4] = 2 * J[2];
        J[5] = 3 * J[1];
        J[6] = 3 * J[2];
    }
    
    {
        F[5] = Vertex::det(a,b) - 0.5 * C0 * ipower(dr0_norm,3);
        array<Real> &J = P[5];
        const Real fac5 = 1.5 * C0 * dr0_norm;
        J[1] =  b.y -  a.x * fac5;
        J[2] = -b.x -  a.y * fac5;
        J[3] = -a.y;
        J[4] =  a.x;
        J[5] = 0;
        J[6] = 0;
    }
    
    {
        F[5] = Vertex::det(a,b) + 3 * Vertex::det(a,c) + 3 * Vertex::det(b,c) - 0.5 * C1 * ipower(dr1_norm,3);
        array<Real> &J = P[6];
        const Real fac6 = 1.5 * C1* dr1_norm;
        J[1] =  3 * c.y +     b.y - fac6 * dr1.x;
        J[2] = -3 * c.x -     b.x - fac6 * dr1.y;
        J[3] =  3 * c.y -     a.y - 2 * fac6 * dr1.x;
        J[4] = -3 * c.x +     a.x - 2 * fac6 * dr1.y;
        J[5] = -3 * b.y - 3 * a.y - 3 * fac6 * dr1.x;
        J[6] =  3 * b.x + 3 * a.x - 3 * fac6 * dr1.y;
    }
}



ArcSolver:: ~ArcSolver() throw() {}

ArcSolver:: ArcSolver() :
P(ARC_NVAR,ARC_NVAR),
U(ARC_NVAR,0.0),
F(ARC_NVAR,0.0),
h(ARC_NVAR,0.0),
ls(ARC_NVAR)
{
    
    
}

#include "yocto/math/kernel/algebra.hpp"
#include "yocto/math/kernel/svd.hpp"

void ArcSolver:: compute(const Arc &arc)
{
    //--------------------------------------------------------------------------
    // initialize
    //--------------------------------------------------------------------------
    arc.init(U);
    
    for(size_t iter=0; iter<30;++iter)
    {
        //-- compute jacobian and value
        arc.estimate(P, F, U);
        std::cerr << "U=" << U << std::endl;
        std::cerr << "P=" << P << std::endl;
        std::cerr << "F=" << F << std::endl;
        
        vector<Real> W(ARC_NVAR,0);
        matrix<Real> V(ARC_NVAR,ARC_NVAR);
        if( !svd<Real>::build(P, W, V) )
        {
            throw exception("Singular SVD arc");
        }
        //std::cerr << "W=" << W << std::endl;
        //std::cerr << "M=" << P << std::endl;
        svd<Real>::truncate(W,1e-3);
        std::cerr << "W=" << W << std::endl;
        //std::cerr << "V=" << V << std::endl;
        svd<Real>::solve(P, W, V, F, h);
        std::cerr << "h=" << h << std::endl;
        
        algebra<Real>::sub(U, h);
        
    }
    
    arc.load(U);
}

