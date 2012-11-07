#include "yocto/utest/run.hpp"
#include "../arc.hpp"



YOCTO_UNIT_TEST_IMPL(arc)
{
    ArcSolver solver;
    Arc arc;
    
    arc.r0 = Vertex(-1,1);
    arc.t0 = Vertex(1,1);
    arc.t0.normalize();
    arc.C0 = -0.3;
    
    arc.r1 = Vertex(1,1);
    arc.t1 = Vertex(1,-1);
    arc.t1.normalize();
    arc.C1 = -0.2;
    
    solver.compute(arc);
    
    ios::ocstream fp("arc.dat",false);
    
    const size_t NP = 100;
    for( size_t i=0; i <= NP; ++i )
    {
        const Real   mu = Real(i)/NP;
        const Vertex v  = arc(mu);
        fp("%g %g\n", v.x, v.y);
    }
    

    
}
YOCTO_UNIT_TEST_DONE()

namespace
{
    class Param
    {
    public:
        inline  Param() {}
        inline ~Param() throw() {}
        
        
        void eval( array<double> &dydx, double x, const array<double> &y )
        {
            
        }
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Param);
    };
}

#include "yocto/math/ode/drvck.hpp"

YOCTO_UNIT_TEST_IMPL(arc2)
{
    Param param;
    ode::field<double>::type eq( &param, & Param::eval );
    ode::drvck<double>::type odeint(1e-4);
    
    odeint.start( 3 );
    
    
}
YOCTO_UNIT_TEST_DONE()
