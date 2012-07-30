#include "yocto/utest/driver.hpp"

using namespace yocto;


YOCTO_UNIT_TEST_INIT(16)
{
    YOCTO_UNIT_TEST_DECL(types);
    YOCTO_UNIT_TEST_DECL(shapes);
    YOCTO_UNIT_TEST_DECL(seg);
    YOCTO_UNIT_TEST_DECL(segpbc);
    YOCTO_UNIT_TEST_DECL(fourier);
    
}
YOCTO_UNIT_TEST_EXEC()
