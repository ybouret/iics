#include "arc.hpp"


Arc:: ~Arc() throw() {}

Arc:: Arc() throw() :
r0(),
t0(),
C0(),
r1(),
t1(),
C1(),
dr()
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


void Arc:: fill( matrix<Real> &P, array<Real> &U ) const
{
    assert( P.cols == 6 );
    assert( P.rows == 6 );
    
    
    dr = r1-r0;
    //--------------------------------------------------------------------------
    // rows 3
    //--------------------------------------------------------------------------
    //P[3][1] = t0.y; P[3][2] = -t0.x;
    
    //--------------------------------------------------------------------------
    // rows 3
    //--------------------------------------------------------------------------
    P[4][1] = t1.y; P[4][2] = -t1.x; P[4][3] = 2*t1.y; P[4][4] = -2*t1.x; P[4][5] = 3*t1.y; P[4][6] = -3*t1.x;
    
#if 1
    const Real len = 0.5*dr.norm();
    a = len * t0;
    U[1] = a.x;
    U[2] = a.y;
    
    const Vertex v1 = dr - a;
    const Vertex v2 = len * t1 - a;
    b = 3.0 * v1 - v2;
    c = v2 - 2.0 * v1;
    
    U[3] = b.x;
    U[4] = b.y;
    
    U[5] = c.x;
    U[6] = c.y;
#else
    
    U[1] = dr.x;
    U[2] = dr.y;
    
    U[3] = 0;
    U[4] = 0;
    U[5] = 0;
    U[6] = 0;
    
#endif
    
}

#include "yocto/code/ipower.hpp"

void Arc:: estimate( matrix<Real> &P, array<Real> &F, const array<Real> &U ) const
{
    load(U);
    const Vertex dot_r0 = a;
    const Vertex dot_r1 = a+2.0*b+3.0*c;
    const Real dot_r0_norm = dot_r0.norm();
    const Real dot_r1_norm = dot_r1.norm();
    
    F[1] = a.x + b.x + c.x - dr.x;        // P[1] = dF[1]/dU was set
    F[2] = a.y + b.y + c.y - dr.y;        // P[2] = dF[2]/dU was set
    
    
    F[3] = dot_r0 * t0 - dot_r0_norm;
    
    {
        // dF[3]/da.x
        P[3][1] = t0.x - a.x / dot_r0_norm;
        
        // dF[3]/da.y
        P[3][2] = t0.y - a.y / dot_r0_norm;
        
        P[3][3] = 0;
        P[3][4] = 0;
        P[3][5] = 0;
        P[3][6] = 0;
        
    }
    
    F[4] = Vertex::det(dot_r1,t1);        // P[4] = dF[4]/dU was set
    
    F[5] = Vertex::det(a,b) - 0.5 * C0 * ipower( dot_r0_norm,3);
    
    {
        const Real fac5 = 1.5 * C0 * dot_r0_norm;
        
        // dF[5]/da.x
        P[5][1] = b.y -  fac5 * a.x;
        
        // dF[5]/da.y
        P[5][2] = -b.x - fac5 * a.y;
        
        // dF[5]/db.x
        P[5][3] = -a.y;
        
        // dF[5]/db.y
        P[5][4] = a.x;
        
        // dF[5]/dc.x
        P[5][5] = 0;
        
        // dF[6]/dc.y
        P[5][6] = 0;
    }
    
    F[6] = Vertex::det(a,b) + 3 * Vertex::det(b,c) + 3 * Vertex::det(a,c) - 0.5 * C1 * ipower( dot_r1_norm, 3 );
    
    {
        const Real fac6 = 1.5 * C1 * dot_r1_norm;
        
        // dF[6]/da.x
        P[6][1] = b.y + 3*c.y    -  fac6 * a.x;
        
        // dF[6]/da.y
        P[6][2] = -b.x - 3*c.x   - fac6 * a.y;
        
        // dF[6]/db.x
        P[6][3] = -a.y + 3*c.y   - fac6 * (2*b.x);
        
        // dF[6]/db.y
        P[6][4] = a.x - 3*c.x    - fac6 * (2*b.y);
        
        // dF[6]/dc.x
        P[6][5] = -3*b.y - 3*a.y - fac6 * (3*c.x);
        
        // dF[6]/dc.y
        P[6][6] = 3*b.x  + 3*a.x - fac6 * (3*c.y);
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
    arc.fill(P,U);
    
    
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

