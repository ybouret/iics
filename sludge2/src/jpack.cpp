#include "jpack.hpp"

JPack:: JPack() throw() :
i(0),
b(0),
y(0),
c(0),
n()
{
}

JPack:: ~JPack() throw() {}

JPack:: JPack( unit_t idx, const Junction *J ) throw() :
i( idx ),
b( J->bubble->id ),
y( J->vertex.y ),
c( J->curvature),
n( J->n )
{
}


JPack:: JPack( const JPack &other ) throw() :
i( other.i ),
b( other.b ),
y( other.y ),
c( other.c ),
n( other.n )
{
}
