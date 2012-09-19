#include "yocto/utest/run.hpp"
#include "../parameters.hpp"


YOCTO_UNIT_TEST_IMPL(param)
{
    const Coord  N(10,20);
    const Vertex L(2.0,3.0);
    Parameters   param(N,L,0,1);
    
}
YOCTO_UNIT_TEST_DONE()

