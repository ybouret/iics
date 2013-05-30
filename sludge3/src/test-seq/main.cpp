#include "yocto/utest/driver.hpp"

using namespace yocto;

YOCTO_UNIT_TEST_INIT(8)
{
    YOCTO_UNIT_TEST_DECL(bubble);
    YOCTO_UNIT_TEST_DECL(grid);
    YOCTO_UNIT_TEST_DECL(curv);
    YOCTO_UNIT_TEST_DECL(segment);

}
YOCTO_UNIT_TEST_EXEC()
