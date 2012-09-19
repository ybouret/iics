#include "yocto/utest/driver.hpp"

using namespace yocto;


YOCTO_UNIT_TEST_INIT(16)
{
    YOCTO_UNIT_TEST_DECL(bubble);
    YOCTO_UNIT_TEST_DECL(bubbles);
    YOCTO_UNIT_TEST_DECL(param);
}
YOCTO_UNIT_TEST_EXEC()
