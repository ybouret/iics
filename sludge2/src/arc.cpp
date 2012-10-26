#include "arc.hpp"


Arc:: ~Arc() throw() {}

Arc:: Arc() throw() :
r0(),
t0(),
C0(),
r1(),
t1(),
C1(),
delta_r()
{
}


Vertex Arc:: operator()( const Real mu ) const throw()
{
    return r0 + mu * a + (mu*mu) * b + (mu*mu*mu) * c;
}

void Arc:: load( const array<Real> &U ) const
{
    assert( U.size() == 6 );
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
    const Vertex dr0 = a;
    const Vertex dr1 = a+2.0*b+3.0*c;
    const Real dr0_norm  = dr0.norm();
    const Real dr0_norm3 = ipower( dr0_norm, 3);
    const Real dr1_norm  = dr1.norm();
    const Real dr1_norm3 = ipower( dr1_norm, 3);

    {
        F[1] = a.x + b.x + c.x - delta_r.x;        // P[1] = dF[1]/dU was set
    }
    
    {
        F[2] = a.y + b.y + c.y - delta_r.y;        // P[2] = dF[2]/dU was set
    }
    
    {
        F[3] = t0.x - dr0.x / dr0_norm;
        
        // dF3/da.x
        P[3][1] = a.x * a.x / dr0_norm3 - 1.0 / dr0_norm;
        
        // dF3/da.y
        P[3][2] = a.x * a.y / dr0_norm3;
        
        P[3][3] = P[3][4] = P[3][5] = P[3][6] = 0;
    }
    
    {
        F[4] = t0.y - dr0.y / dr0_norm;
        
        // dF4/da.x
        P[4][1] = a.x * a.y / dr0_norm3;
        
        // dF4/da.y
        P[4][2] = a.y*a.y / dr0_norm3 - 1.0 / dr0_norm;
        
        P[4][3] = P[4][4] = P[4][5] = P[4][6] = 0;
    }
    
    {
        F[5] = t1.x - dr1.x / dr1_norm;
        
        P[5][1] = dr1.x * dr1.x / dr1_norm3 - 1.0/dr1_norm;
        P[5][2] = dr1.x * dr1.y / dr1_norm3;
        P[5][3] = 2 * P[5][1];
        P[5][4] = 2 * P[5][2];
        P[5][5] = 3 * P[5][1];
        P[5][6] = 3 * P[5][2];
    }
    
    {
        F[6] = t1.y - dr1.y / dr1_norm;
        P[6][1] = P[5][2];
        P[6][2] = dr1.y * dr1.y / dr1_norm3 - 1.0/dr1_norm;
        P[6][3] = 2 * P[6][1];
        P[6][4] = 2 * P[6][2];
        P[6][5] = 3 * P[6][1];
        P[6][6] = 3 * P[6][2];
    }
}



ArcSolver:: ~ArcSolver() throw() {}

ArcSolver:: ArcSolver() :
P(6,6),
iP(6,6),
U(6,0.0),
F(6,0.0),
h(6,0.0),
ls(6)
{
    //--------------------------------------------------------------------------
    // first two rows: won't change
    //--------------------------------------------------------------------------
    P[1][1] = 1; P[1][2] = 0; P[1][3] = 1; P[1][4] = 0; P[1][5] = 1; P[1][6] = 0;
    P[2][1] = 0; P[2][2] = 1; P[2][3] = 0; P[2][4] = 1; P[2][5] = 0; P[2][6] = 1;
    
}

#include "yocto/math/kernel/algebra.hpp"

void ArcSolver:: compute(const Arc &arc)
{
    //--------------------------------------------------------------------------
    // initialize
    //--------------------------------------------------------------------------
    arc.init(U);
        
    for(size_t iter=0; iter<10;++iter)
    {
        //-- compute jacobian and value
        arc.estimate(P, F, U);
        std::cerr << "U=" << U << std::endl;
        std::cerr << "P=" << P << std::endl;
        std::cerr << "F=" << F << std::endl;
        
        //-- compute Newton's step
        iP.assign(P);
        if( !ls.LU(iP) )
        {
            throw exception("singular arc");
        }
        for( size_t i=1; i <= U.size(); ++i ) h[i] = -F[i];
        ls(iP,h);
        std::cerr << "h=" << h << std::endl;
        
        algebra<Real>::add(U, h);
        
    }
    
    arc.load(U);
}

