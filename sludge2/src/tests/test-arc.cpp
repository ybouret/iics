#include "yocto/utest/run.hpp"
#include "../arc.hpp"


YOCTO_UNIT_TEST_IMPL(arc)
{
    ArcSolver solver;
    
    Arc arc;
    
    arc.r0 = Vertex(-1,1);
    arc.t0 = Vertex(1,1);
    arc.t0.normalize();
    arc.C0 = -1;
    
    
    arc.r1 = Vertex(1,1);
    arc.t1 = Vertex(1,-0.5);
    arc.t1.normalize();
    arc.C1 = 0;
    
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
