#include "yocto/utest/run.hpp"
#include "../bubbles.hpp"

YOCTO_UNIT_TEST_IMPL(bubble)
{
    const PBC pbc(1);
    Bubbles   bubbles(pbc);
    bubbles.create(10);
    
}
YOCTO_UNIT_TEST_DONE()
