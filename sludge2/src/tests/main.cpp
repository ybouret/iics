#include "yocto/utest/driver.hpp"

using namespace yocto;


YOCTO_UNIT_TEST_INIT(16)
{
	YOCTO_UNIT_TEST_DECL(types);
	YOCTO_UNIT_TEST_DECL(tracers);
    YOCTO_UNIT_TEST_DECL(spots);
    YOCTO_UNIT_TEST_DECL(bubble);
    YOCTO_UNIT_TEST_DECL(seg);
    YOCTO_UNIT_TEST_DECL(param);
	YOCTO_UNIT_TEST_DECL(arc);
}
YOCTO_UNIT_TEST_EXEC()
