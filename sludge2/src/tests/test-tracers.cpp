#include "yocto/utest/run.hpp"
#include "../tracer.hpp"

YOCTO_UNIT_TEST_IMPL(tracers)
{
    Tracer::Cache tcache;
    Tracers       tracers( tcache );
    
    std::cerr << "sizeof(Tracer)=" << sizeof(Tracer) << std::endl;
    
    for( size_t i=0; i <10; ++i ) tracers.append();
    std::cerr << "#tracers=" << tracers.size << std::endl;
    tracers.empty();
    for( size_t i=0; i <15; ++i ) tracers.append();
    std::cerr << "#tracers=" << tracers.size << std::endl;

}
YOCTO_UNIT_TEST_DONE()

