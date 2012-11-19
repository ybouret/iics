#include "arc.hpp"

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

Arc:: Arc(const ArcPoint &a_A,
          const ArcPoint &a_Q,
          const ArcPoint &a_B ) throw() :
A( a_A ),
Q( a_Q ),
B( a_B ),
AQ( A.r, Q.r ),
QB( Q.r, B.r ),
AB( A.r, B.r ),
lambda( clamp<Real>(0,(AQ*AB)/(AB*AB),1) ) ,
omlam(1-lambda),
lambda2(lambda*lambda),
lambda3(lambda2*lambda)
{
    A.delta = Vertex::angle_of( A.t, Q.t);
    Q.delta = Vertex::angle_of( Q.t, B.t);
    std::cerr << "lambda=" << lambda << std::endl;
    //std::cerr << "A.theta=" << A.theta << " -> Q.theta=" << Q.theta <<  ": " << A.delta << std::endl;
    //std::cerr << "Q.theta=" << Q.theta << " -> B.theta=" << B.theta <<  ": " << Q.delta << std::endl;
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
            return Q.theta;
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
            return B.theta;
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
    const Real ftol = 1e-5;
    K.ldz();
    K[1][1] = 1; K[1][3] = 1; K[1][5] = 1;
    K[2][2] = 1; K[2][4] = 1; K[2][6] = 1;
    K[3][1] = lambda; K[3][3] = lambda2; K[3][5] = lambda3;
    
    U[1] = B.alpha - A.alpha;
    U[2] = B.beta  - A.beta;
    U[3] = Q.alpha - A.alpha;
    
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
    
    {
        numeric<Real>::function dI( this, & Arc::dI_QB);
        numeric<Real>::function dJ( this, & Arc::dJ_QB);
        nu = 0;
        const Real I0 = integrate<Real>(lambda,1, dI, ftol);
        const Real J0 = integrate<Real>(lambda,1, dJ, ftol);
        U[5] = (B.P - Q.P) - (I0*Q.alpha+J0*Q.beta);
        
        nu=1;
        K[5][1] = integrate<Real>(lambda,1, dI, ftol);
        K[5][2] = integrate<Real>(lambda,1, dJ, ftol);
        
        nu=2;
        K[5][3] = integrate<Real>(lambda,1, dI, ftol);
        K[5][4] = integrate<Real>(lambda,1, dJ, ftol);
        
        nu=3;
        K[5][5] = integrate<Real>(lambda,1, dI, ftol);
        K[5][6] = integrate<Real>(lambda,1, dJ, ftol);
    }

    
    
    
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
J(NX,NX),
JK(NX,NC),
H(NC,NC),
U(NC,0),
W(NC,0),
V(NC,NC),
Lambda(NC,0),
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
    arc.load(K,U);
       
    
    std::cerr << "K=" << K << std::endl;
    mkl::mul_rtrn(JK, J,K);
    std::cerr << "JK=" << JK << std::endl;
    mkl::mul(H,K,JK);
    std::cerr << "H=" << H << std::endl;
    std::cerr << "U="<< U << std::endl;
    if( !svd<Real>::build(H, W, V))
        throw exception("Invalid Arc");
    std::cerr << "W=" << W << std::endl;
    svd<Real>::truncate(W, 1e-9);
    svd<Real>::solve(H, W, V, U, Lambda);
    std::cerr << "Lambda=" << Lambda << std::endl;
    mkl::mul(X,JK,Lambda);
    std::cerr << "X=" << X << std::endl;
    const Real gQ = arc.A.beta + arc.lambda * X[4] + arc.lambda2 * X[5] + arc.lambda3 * X[6];
    std::cerr << "gQ=" << gQ << std::endl;
}



