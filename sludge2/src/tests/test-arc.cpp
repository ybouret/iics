#include "yocto/utest/run.hpp"
#include "../arc.hpp"

#include "yocto/math/kernel/svd.hpp"
#include "yocto/math/kernel/algebra.hpp"

#include "yocto/string/conv.hpp"

YOCTO_UNIT_TEST_IMPL(arc)
{
    
    Real lambda = 0.5;
    if( argc >1 )
        lambda = strconv::to_real<Real>( argv[1], "lambda" );
    matrix<Real> H(4,4);
    vector<Real> W(4,0);
    matrix<Real> V(4,4);
    vector<Real> L(4,0);
    vector<Real> U(4,0);
    matrix<Real> JK(4,4);
    vector<Real> X(4,0);
    
    const Vertex A( 10*(0.5-Alea()), 10*(0.5-Alea()) );
    const Vertex B( 10*(0.5-Alea()), 10*(0.5-Alea()) );
    const Vertex AB(A,B);
    const Vertex AQ( lambda*AB.x, lambda*(1-lambda)*10*(0.5-Alea()));
    
    const Real lambda2 = lambda*lambda;
    const Real lambda3 = lambda2 * lambda;
    const Real lambda4 = lambda2 * lambda2;
    const Real lfac    = 3.0 * lambda4 - 6.0*lambda3 + 4.0 * lambda2;
    
    H.ldz();
    H[1][1] = 1; H[1][3] = lambda;
    H[2][2] = 1; H[2][4] = lambda;
    H[3][1] = lambda; H[3][3] = lfac;
    H[4][2] = lambda; H[4][4] = lfac;
    
    JK.ldz();
    JK[1][1] = 1.0;
    JK[2][2] = 1.0;
    JK[1][3] = JK[2][4] = 4*lambda -3*lambda2;
    JK[3][3] = JK[4][4] = 3*(lambda2-lambda);
    std::cerr << "JK=" << JK << std::endl;
    std::cerr << "H=" << H << std::endl;
    if( !svd<Real>::build( H, W, V ) )
    {
        throw exception("Invalid Arc");
    }
    std::cerr << "W0=" << W << std::endl;
    (void) svd<Real>::truncate(W, numeric<double>::ftol);
    std::cerr << "W1=" << W << std::endl;
    
    U[1] = AB.x;
    U[2] = AB.y;
    U[3] = AQ.x;
    U[4] = AQ.y;
    //for(size_t i=4;i>0;--i) U[i] *= 6;
    
    svd<Real>::solve(H,W,V,U,L);
    algebra<Real>::mul(X, JK, L);
    //for(size_t i=4;i>0;--i) X[i]/=6;
    std::cerr << "X=" << X << std::endl;
    std::cerr << "AB=" << AB << std::endl;
    const Vertex a(X[1],X[2]);
    const Vertex b(X[3],X[4]);
    std::cerr << "a=" << a << std::endl;
    std::cerr << "b=" << b << std::endl;
    std::cerr << "a+b=" << a+b << std::endl;
    
    std::cerr << "AQ=" << AQ << std::endl;
    std::cerr << "lambda*a+lambda^2*b=" << lambda *a + lambda2 *b << std::endl;
    
    
    {
        ios::ocstream fp("arc.dat",false);
        for(size_t i=0;i<=100;++i)
        {
            const Real   mu = Real(i)/100;
            const Vertex r  = A + mu * a + mu*mu*b;
            fp("%g %g\n",r.x,r.y);
        }
        
    }
    
    {
        ios::ocstream fp("tri.dat",false);
        const Vertex Q = A + AQ;
        fp("%g %g\n", A.x, A.y);
        fp("%g %g\n", Q.x, Q.y);
        fp("%g %g\n", B.x, B.y);
    }
#if 0
    ArcSolver solver;
    
    {
        // Circle with same pressure
        ArcPoint A( Vertex(-2,0), Vertex(0,1),  0.2, 0, 1 );
        ArcPoint Q( Vertex(0,2),  Vertex(1,0),  0.2, 0, 0 );
        ArcPoint B( Vertex(2,0),  Vertex(0,-1), 0.2, 0, 1 );
        
        Arc arc(A,Q,B);
        
        
        solver(arc);
        std::cerr << "gA=" << A.beta << ", gB=" << B.beta << " => gQ=" << Q.beta << std::endl;
    }
    
    
    {
        // Pressure ramp: beta=1, alpha=0
        ArcPoint A( Vertex(0,0),   Vertex(1,0), 0.2, 0, 1);
        ArcPoint Q( Vertex(1.2,0), Vertex(1,0), 1.4, 0, 0);
        ArcPoint B( Vertex(2.0,0), Vertex(1,0), 2.2, 0, 1);
        
        Arc arc(A,Q,B);
        solver(arc);
        std::cerr << "gA=" << A.beta << ", gB=" << B.beta << " => gQ=" << Q.beta << std::endl;
    }
#endif
    
    
}
YOCTO_UNIT_TEST_DONE()
