#include "yocto/utest/run.hpp"
#include "../tracer.hpp"

YOCTO_UNIT_TEST_IMPL(tracers)
{
    Tracer::Cache tcache;
    Tracers       tracers( tcache );
    for( size_t i=0; i <10; ++i ) tracers.append();
    std::cerr << "#tracers=" << tracers.size << std::endl;
    
}
YOCTO_UNIT_TEST_DONE()

