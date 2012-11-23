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
delta(0)
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
#include "yocto/code/utils.hpp"

#include "yocto/math/kernel/algebra.hpp"
#include "yocto/code/ipower.hpp"
////////////////////////////////////////////////////////////////////////////////
//
// Arc
//
////////////////////////////////////////////////////////////////////////////////

Arc:: ~Arc() throw() {}

static inline
Real compute_lambda( const Vertex &AQ, const Vertex &AB )
{
    const Real num = AQ*AB;
    if( num <= 0 )
        return 0;
    else
    {
        const Real den = AB*AB;
        if( num >= den )
            return 1;
        else
            return num/den;
    }
}

Arc:: Arc(const ArcPoint &a_A,
          const ArcPoint &a_Q,
          const ArcPoint &a_B ) throw() :
A( a_A ),
Q( a_Q ),
B( a_B ),
AQ( A.r, Q.r ),
QB( Q.r, B.r ),
AB( A.r, B.r ),
lambda( compute_lambda(AQ, AB) ) ,
omlam(1-lambda),
lambda2(lambda*lambda),
lambda3(lambda2*lambda)
{
    A.delta = Vertex::angle_of( A.t, Q.t);
    Q.delta = Vertex::angle_of( Q.t, B.t);
    
    //std::cerr << "A=" << A.r << ", Q=" << Q.r << ", B=" << B.r << "AQ=" << AQ << ", AB=" << AB << ", lambda=" << lambda << std::endl;
    std::cerr << "#A.theta=" << A.theta << " -> Q.theta=" << Q.theta <<  ": " << A.delta << std::endl;
    std::cerr << "#Q.theta=" << Q.theta << " -> B.theta=" << B.theta <<  ": " << Q.delta << std::endl;
}

Real Arc:: theta_AQ(Real mu) const throw()
{
    if(mu<=0)
    {
        return A.theta;
    }
    else
    {
        if(mu>=lambda)
        {
            return A.theta+A.delta;
        }
        else
        {
            return A.theta + (mu * A.delta) / lambda;
        }
    }
}

Real Arc:: theta_QB(Real mu) const throw()
{
    if(mu<=lambda)
    {
        return Q.theta;
    }
    else
    {
        if(mu>=1)
        {
            return Q.theta + Q.delta;
        }
        else
        {
            return Q.theta + ( (mu-lambda) * Q.delta ) / omlam;
        }
    }
}

Real Arc:: theta(Real mu) const throw()
{
    return mu <= lambda ? theta_AQ(mu) : theta_QB(mu);
}


Vertex Arc:: rdot_AQ(Real) const throw()
{
    return AQ/lambda;
}

Vertex Arc:: rdot_QB(Real) const throw()
{
    return QB/omlam;
}

void Arc:: load(matrix<Real> &K, array<Real> &U) const throw()
{
    const Real ftol = 1e-7;
    K.ldz();
    K[1][1] = 1; K[1][3] = 1; K[1][5] = 1;
    K[2][2] = 1; K[2][4] = 1; K[2][6] = 1;
    K[3][1] = lambda; K[3][3] = lambda2; K[3][5] = lambda3;
    
    U[1] = B.alpha - A.alpha;
    U[2] = B.beta  - A.beta;
    U[3] = Q.alpha - A.alpha;
    
    if(lambda>0)
    {
        numeric<Real>::function dI( this, & Arc::dI_AQ);
        numeric<Real>::function dJ( this, & Arc::dJ_AQ);
        nu = 0;
        const Real I0 = integrate<Real>(0,lambda, dI, ftol);
        const Real J0 = integrate<Real>(0,lambda, dJ, ftol);

        U[4] = (Q.P - A.P) - (I0*A.alpha+J0*A.beta);
        
        nu=1;
        K[4][1] = integrate<Real>(0,lambda, dI, ftol);
        K[4][2] = integrate<Real>(0,lambda, dJ, ftol);
        
        nu=2;
        K[4][3] = integrate<Real>(0,lambda, dI, ftol);
        K[4][4] = integrate<Real>(0,lambda, dJ, ftol);
        
        nu=3;
        K[4][5] = integrate<Real>(0,lambda, dI, ftol);
        K[4][6] = integrate<Real>(0,lambda, dJ, ftol);
    }
    else U[4] = 0;
    
    if(lambda<1)
    {
        numeric<Real>::function dI( this, & Arc::dI_QB);
        numeric<Real>::function dJ( this, & Arc::dJ_QB);
        nu = 0;
        const Real I0 = integrate<Real>(1,lambda, dI, ftol);
        const Real J0 = integrate<Real>(1,lambda, dJ, ftol);
       
        U[5] = (Q.P - B.P) - (I0*B.alpha+J0*B.beta);
        
        nu=1;
        K[5][1] = integrate<Real>(1,lambda,dI, ftol);
        K[5][2] = integrate<Real>(1,lambda,dJ, ftol);
        
        nu=2;
        K[5][3] = integrate<Real>(1,lambda,dI, ftol);
        K[5][4] = integrate<Real>(1,lambda,dJ, ftol);
        
        nu=3;
        K[5][5] = integrate<Real>(1,lambda,dI, ftol);
        K[5][6] = integrate<Real>(1,lambda,dJ, ftol);
    }
    else
        U[5] = 0;
    
    
    
    
}


Real Arc:: dI_AQ(Real mu) const throw()
{
    const Real   th = theta_AQ(mu);
    const Vertex v  = rdot_AQ(mu);
    const Vertex tau( Cos(th), Sin(th) );
    return ipower(mu,nu) * (v*tau);
}

Real Arc:: dJ_AQ(Real mu) const throw()
{
    const Real   th = theta_AQ(mu);
    const Vertex v  = rdot_AQ(mu);
    const Vertex n( -Sin(th), Cos(th) );
    return ipower(mu,nu) * (v*n);
}

Real Arc:: dI_QB(Real mu) const throw()
{
    const Real   th = theta_QB(mu);
    const Vertex v  = rdot_QB(mu);
    const Vertex tau( Cos(th), Sin(th) );
    return ipower(mu,nu) * (v*tau);
}

Real Arc:: dJ_QB(Real mu) const throw()
{
    const Real   th = theta_QB(mu);
    const Vertex v  = rdot_QB(mu);
    const Vertex n( -Sin(th), Cos(th) );
    return ipower(mu,nu) * (v*n);
}


////////////////////////////////////////////////////////////////////////////////
//
// ArcSolver
//
////////////////////////////////////////////////////////////////////////////////
ArcSolver:: ~ArcSolver() throw() {}

ArcSolver:: ArcSolver() :
K(NC,NX),
U(NX,NC),
W(NC,0),
V(NC,NC),
Z(NC,0),
J(NX,NX),
JU(NX,NC),
H(NC,NC),
Y(NC,0),
solve(NC),
X(NX,0)
{
    static const Real __J[] =
    {
        9,0,-18,0,10,0,
        0,9,0,-18,0,10,
        -18,0,48,0,-30,0,
        0,-18,0,48,0,-30,
        10,0,-30,0,20,0,
        0,10,0,-30,0,20
    };
    assert(sizeof(__J)==J.cols*J.rows*sizeof(Real));
    memcpy( &J[1][1], __J, sizeof(__J) );
    std::cerr << "J=" << J << std::endl;
}


typedef algebra<Real> mkl;

void ArcSolver:: operator()( const Arc &arc )
{
    //-- get the matrix/vector
    arc.load(K,Z);
    //std::cerr << "K=" << K << std::endl;
    //std::cerr << "Z=" << Z << std::endl;
    
    //-- build the svd
    for( size_t i=1; i <= NC; ++i )
        for(size_t j=1; j <= NX; ++j ) U[j][i] = K[i][j];
    if( ! svd<Real>::build(U, W, V) )
        throw exception("Invalid Arc, level-1");
    //std::cerr << "U=" << U << std::endl;
    //std::cerr << "W=" << W << std::endl;
    //std::cerr << "V=" << V << std::endl;
    
    //-- JU= J*U
    mkl::mul(JU, J, U);
    //std::cerr << "JU=" << JU << std::endl;
    
    //-- H = U'*J*U
    mkl::mul_ltrn(H, U, JU);
    //std::cerr << "H=" << H << std::endl;
    if( ! solve.LU(H) )
        throw exception("Invalid Arc, level-2");
    
    svd<Real>::truncate(W, 1e-8);
    matrix<Real> W0; W0.diag(W);
    matrix<Real> iW(NC,NC);
    for( size_t i=1; i <= NC; ++i)
        if( Fabs(W[i])>0 ) iW[i][i] = 1.0/W[i];
    //std::cerr << "W0=" << W0 << std::endl;
    //std::cerr << "iW=" << iW << std::endl;
    
    //-- Y1 = inv(W)*V'*Z
    mkl::mul_trn(Y, V, Z);
    for( size_t i=NC;i>0;--i)
    {
        const Real Wi = W[i];
        if( Fabs(Wi)>0 )
            Y[i] /= Wi;
        else
            Y[i] = 0;
    }
    //std::cerr << "Y1=" << Y << std::endl;
    
    //-- Y2 = inv(H) * Y1
    solve(H,Y);
    //std::cerr << "Y2=" << Y << std::endl;

    //-- Y3 = W*inv(W) * Y2
    for( size_t i=NC;i>0;--i)
    {
        if( Fabs(W[i]) <= 0 )
            Y[i] = 0;
    }
    //std::cerr << "Y3=" << Y << std::endl;

    mkl::mul(X, JU, Y);
    //std::cerr << "X=" << X << std::endl;

    const Real gQ = arc.A.beta + arc.lambda * X[4] + arc.lambda2 * X[5] + arc.lambda3 * X[6];
    //std::cerr << "gQ=" <<  gQ << std::endl;
    (Real &)(arc.Q.beta) = gQ;
}
#endif



