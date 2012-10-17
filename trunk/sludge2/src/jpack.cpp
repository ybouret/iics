#include "jpack.hpp"

JPack:: JPack() throw() :
i(0),
b(0),
y(0),
c(0),
t(),
n(),
g(0)
{
}

JPack:: ~JPack() throw() {}

JPack:: JPack( unit_t idx, const Junction *J ) throw() :
i( idx ),
b( J->bubble->id ),
y( J->vertex.y ),
c( J->curvature),
t( J->t ),
n( J->n ),
g( J->gradP_t)
{
}


JPack:: JPack( const JPack &other ) throw() :
i( other.i ),
b( other.b ),
y( other.y ),
c( other.c ),
t( other.t ),
n( other.n ),
g( other.g )
{
}
