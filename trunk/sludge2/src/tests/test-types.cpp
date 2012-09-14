#include "yocto/utest/run.hpp"
#include "../types.hpp"

YOCTO_UNIT_TEST_IMPL(types)
{
    PBC pbc(2.2);
    for( Real y=-5;y<=5; y += 0.1)
    {
        std::cerr << "y=" << y << " => " << pbc.apply(y) << std::endl;
    }
    
}
YOCTO_UNIT_TEST_DONE()
