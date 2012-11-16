#include "yocto/utest/run.hpp"
#include "../arc.hpp"

YOCTO_UNIT_TEST_IMPL(arc)
{
    ArcSolver solver;
    
    {
        // Circle with same pressure
        ArcPoint A( Vertex(-2,0), Vertex(0,1),  0.2, 0, 1 );
        ArcPoint Q( Vertex(0,2),  Vertex(1,0),  0.2, 0, 0 );
        ArcPoint B( Vertex(2,0),  Vertex(0,-1), 0.2, 0, 1 );
        
        Arc arc(A,Q,B);
        
        
        solver(arc);
    }
    
    
    {
        // Pressure ramp: beta=1, alpha=0
        ArcPoint A( Vertex(0,0),   Vertex(1,0), 0.2, 0, 1);
        ArcPoint Q( Vertex(1.2,0), Vertex(1,0), 1.4, 0, 0);
        ArcPoint B( Vertex(2.0,0), Vertex(1,0), 2.2, 0, 1);
        
        Arc arc(A,Q,B);
        solver(arc);
    }
    
    
    
}
YOCTO_UNIT_TEST_DONE()
