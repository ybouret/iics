#include "yocto/utest/run.hpp"
#include "../spot.hpp"

YOCTO_UNIT_TEST_IMPL(spots)
{
    Tracer::Cache tcache;
    Tracers       tracers( tcache );
    Spot::Cache   scache;
    Spots         spots( scache );
    
    for( size_t i=0; i < 20; ++i ) tracers.append();
    
    for( size_t i=0; i <10; ++i)
        spots.attach( tracers.fetch(i) );
    

}
YOCTO_UNIT_TEST_DONE()
