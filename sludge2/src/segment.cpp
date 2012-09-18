#include "segment.hpp"

Segment:: ~Segment() throw() {}

Segment:: Segment( const Real x_or_y, Junction::Cache &jcache ) throw() :
Junctions( jcache ),
value( x_or_y)
{
}
