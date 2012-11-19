#include "arc.hpp"

#if 0
////////////////////////////////////////////////////////////////////////////////
//
// ArcPoint
//
////////////////////////////////////////////////////////////////////////////////

ArcPoint:: ArcPoint(const Vertex &a_r,
                    const Vertex &a_t,
                    const Real    a_P,
                    const Real    gt,
                    const Real    gn) :
r( a_r ),
t( a_t ),
P( a_P ),
alpha(gt),
beta(gn),
theta( t.angle() ),
delta(0),
C0(0),
S0(0),
C1(0),
S1(0),
C2(0),
S2(0)
{
    ((Vertex &)t).normalize();
}

#if 0
ArcPoint:: ArcPoint( const ArcPoint &other ) throw() :
r( other.r ),
t( other.t ),
P( other.P ),
alpha( other.alpha ),
beta(  other.beta  ),
theta( other.theta ),
delta( other.delta )
{
    
}
#endif

ArcPoint:: ~ArcPoint() throw() {}

#include "yocto/math/fcn/integrate.hpp"

Real ArcPoint:: dC0( Real mu ) const
{
    return Cos( theta + mu * delta );
}

Real ArcPoint:: dS0( Real mu ) const
{
    return Sin( theta + mu * delta );
}

Real ArcPoint:: dC1( Real mu ) const
{
    return mu*Cos( theta + mu * delta );
}

Real ArcPoint:: dS1( Real mu ) const
{
    return mu*Sin( theta + mu * delta );
}

Real ArcPoint:: dC2( Real mu ) const
{
    return mu*mu*Cos( theta + mu * delta );
}

Real ArcPoint:: dS2( Real mu ) const
{
    return mu*mu*Sin( theta + mu * delta );
}

void ArcPoint::computeIntegrals() const
{
    const Real ftol = 1e-5;
    {
        numeric<Real>::function F( this, & ArcPoint::dC0 );
        C0 = integrate<Real>(0,1,F,ftol);
        //std::cerr << "\tC0=" << C0 << std::endl;
    }
    
    {
        numeric<Real>::function F( this, & ArcPoint::dS0 );
        S0 = integrate<Real>(0,1,F,ftol);
        //std::cerr << "\tS0=" << S0 << std::endl;
    }
    
    {
        numeric<Real>::function F( this, & ArcPoint::dC1 );
        C1 = integrate<Real>(0,1,F,ftol);
        //std::cerr << "\tC1=" << C1 << std::endl;
    }
    
    {
        numeric<Real>::function F( this, & ArcPoint::dS1 );
        S1 = integrate<Real>(0,1,F,ftol);
        //std::cerr << "\tS1=" << S1 << std::endl;
    }
    
    {
        numeric<Real>::function F( this, & ArcPoint::dC2 );
        C2 = integrate<Real>(0,1,F,ftol);
       // std::cerr << "\tC2=" << C2 << std::endl;
    }
    
    {
        numeric<Real>::function F( this, & ArcPoint::dS2 );
        S2 = integrate<Real>(0,1,F,ftol);
        //std::cerr << "\tS2=" << S2 << std::endl;
    }
    
}

////////////////////////////////////////////////////////////////////////////////
//
// Arc
//
////////////////////////////////////////////////////////////////////////////////

Arc:: ~Arc() throw() {}

Arc:: Arc(const ArcPoint &a_A,
          const ArcPoint &a_Q,
          const ArcPoint &a_B ) throw() :
A( a_A ),
Q( a_Q ),
B( a_B ),
AQ( A.r, Q.r ),
BQ( B.r, Q.r ),
I0(0),
J0(0),
I1(0),
J1(0),
I2(0),
J2(0)
{
    //--------------------------------------------------------------------------
    // Finalize Data
    //--------------------------------------------------------------------------
    A.delta = Vertex::angle_of( A.t, Q.t );
    B.delta = Vertex::angle_of( B.t, Q.t );
    
    //--------------------------------------------------------------------------
    // Compute Integrals
    //--------------------------------------------------------------------------
    //std::cerr << "A:" << std::endl;
    //std::cerr << "\ttheta=" << A.theta << std::endl;
    //std::cerr << "\tdelta=" << A.delta << std::endl;
    
    //std::cerr << "Compute A Integrals" << std::endl;
    A.computeIntegrals();
    
    //std::cerr << "B:" << std::endl;
    //std::cerr << "\ttheta=" << B.theta << std::endl;
    //std::cerr << "\tdelta=" << B.delta << std::endl;
    
    //std::cerr << "Compute B Integrals" << std::endl;
    B.computeIntegrals();
    
    //--------------------------------------------------------------------------
    // Compute pressure integrals
    //--------------------------------------------------------------------------
    //std::cerr << "Compute Pressure Integrals" << std::endl;
    
    // AQ
    I0 =  A.C0 * AQ.x + A.S0 * AQ.y;
    J0 = -A.S0 * AQ.x + A.C0 * AQ.y;
    
    I1 =  A.C1 * AQ.x + A.S1 * AQ.y;
    J1 = -A.S1 * AQ.x + A.C1 * AQ.y;
    
    I2 =  A.C2 * AQ.x + A.S2 * AQ.y;
    J2 = -A.S2 * AQ.x + A.C2 * AQ.y;
    
    Real tmp = 0;
    
    tmp = I1 - I2;
    eta = I1*I1 + tmp*tmp;
    tmp = J1 - J2;
    eta += J1*J1 + tmp*tmp;
    
    
    
    // BQ
    I0p =  B.C0 * BQ.x + B.S0 * BQ.y;
    J0p = -B.S0 * BQ.x + B.C0 * BQ.y;
    
    I1p =  B.C1 * BQ.x + B.S1 * BQ.y;
    J1p = -B.S1 * BQ.x + B.C1 * BQ.y;
    
    I2p =  B.C2 * BQ.x + B.S2 * BQ.y;
    J2p = -B.S2 * BQ.x + B.C2 * BQ.y;
    
    tmp   = I1p - I2p;
    etap  = I1p*I1p + tmp*tmp;
    tmp   = J1p - J2p;
    etap += J1p*J1p + tmp*tmp;
    
#if 0
    std::cerr << "I0   = " << I0 << std::endl;
    std::cerr << "J0   = " << J0 << std::endl;
    std::cerr << "I1   = " << I1 << std::endl;
    std::cerr << "J1   = " << J1 << std::endl;
    std::cerr << "I2   = " << I2 << std::endl;
    std::cerr << "J2   = " << J2 << std::endl;
    
    std::cerr << "I0p  = " << I0p << std::endl;
    std::cerr << "J0p  = " << J0p << std::endl;
    std::cerr << "I1p  = " << I1p << std::endl;
    std::cerr << "J1p  = " << J1p << std::endl;
    std::cerr << "I2p  = " << I2p << std::endl;
    std::cerr << "J2p  = " << J2p << std::endl;
    
    std::cerr << "eta  =" << eta << std::endl;
    std::cerr << "etap =" << etap << std::endl;
#endif
    
}

void Arc:: load( matrix<Real> &H, array<Real> &U , matrix<Real> &JK) const throw()
{
    //-- build H * 6
    H.ldz();
    H[1][1] = 1;   H[1][4] =  I1;
    H[2][2] = 1;   H[2][5] =  I1p;
    H[3][3] = 2;   H[3][4] =  J1;  H[3][5] = -J1p;
    H[4][1] = I1;  H[4][3] =  J1;  H[4][4] = eta;
    H[5][2] = I1p; H[5][3] = -J1p; H[5][5] = etap;
    //std::cerr << "H=" << H << "/6" << std::endl;
   // std::cerr << "H6=" << H  << std::endl;

    //-- build unknown vector
    U[1] = Q.alpha - A.alpha;
    U[2] = Q.alpha - B.alpha;
    U[3] = B.beta  - A.beta;
    U[4] = (Q.P - A.P) - ( A.alpha * I0  + A.beta * J0 );
    U[5] = (Q.P - B.P) - ( B.alpha * I0p + B.beta * J0p);
    //std::cerr << "U=" << U << std::endl;
    
    for( size_t i=U.size();i>0;--i)
        U[i] *= 6;
    
    JK.ldz();
    JK[1][1] = JK[5][2] = JK[3][3] = 1.0/6;
    JK[7][3] = -1.0/6;
    
    
    JK[1][4] = (2*I1)/3 - I2/2;
    JK[2][4] = (I2-I1)/2;
    JK[3][4] = (2*J1)/3 - J2/2;
    JK[4][4] = (J2-J1)/2;
    
    JK[5][5] = (2*I1p)/3 - I2p/2;
    JK[6][5] = (I2p-I1p)/2;
    JK[7][5] = (2*J1p)/3 - J2p/2;
    JK[8][5] = (J2p-J1p)/2;
    
    //std::cerr << "JK=" << JK << std::endl;
    
    
    
}


////////////////////////////////////////////////////////////////////////////////
//
// ArcSolver
//
////////////////////////////////////////////////////////////////////////////////
#include "yocto/math/kernel/algebra.hpp"
ArcSolver:: ArcSolver() :
H(N,N),
U(N,0),
L(N,0),
W(N,0),
V(N,N),
JK(M,N),
X(M,0)
{
}

ArcSolver:: ~ArcSolver() throw()
{
}

void ArcSolver:: operator()( const Arc &arc )
{
    //--------------------------------------------------------------------------
    // create the matrix and the vector to compute lagrange multipliers
    //--------------------------------------------------------------------------
    arc.load(H,U,JK);

    if( ! svd<Real>::build(H,W,V) )
        throw exception("Invalid Arc");
    //std::cerr << "W=" << W << std::endl;
    //--------------------------------------------------------------------------
    // TODO: truncate
    // compute the lagrange multipliers into L
    //--------------------------------------------------------------------------
    svd<Real>::solve(H, W, V, U, L);
    //std::cerr << "L=" << L << std::endl;
    
    //--------------------------------------------------------------------------
    // compute the solution vector
    //--------------------------------------------------------------------------
    algebra<double>::mul(X, JK, L);
    
    //std::cerr << "X=" << X << std::endl;
    const Real gQ1 = arc.A.beta + X[3] + X[4];
    const Real gQ2 = arc.B.beta + X[7] + X[8];

    std::cerr << "gQ1=" << gQ1 << std::endl;
    std::cerr << "gQ2=" << gQ2 << std::endl;
    
    (Real &)(arc.Q.beta) = (gQ1+gQ2)*0.5;
}
#endif


