#include "junctions.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

Junction:: Junction( List &r) throw() :
prev(0),
next(0),
root(r)
{
}

Junction:: ~Junction() throw()
{
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Junction::List:: List( Type t, const Real &v) throw() :
type(t),
level(v)
{
}

#if 0
Junction::List:: List(  const List &other ) throw() :
type( other.type ),
level( other.level )
{
}
#endif

Junction::List:: ~List() throw() { auto_delete(); }

Junction * Junction::List::append()
{
    Junction *J = new Junction(*this);
    push_back(J);
    return J;
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
#include "yocto/memory/global.hpp"

Junctions:: Junctions( const Grid &grid ) :
Layout( grid ),
jcount( grid.width.x + grid.width.y ),
jlists( memory::kind<memory::global>::acquire_as<Junction::List>(jcount) ),
jvert(jlists),
jhorz(jvert+width.x)
{
    jvert -= lower.x;
    jhorz -= lower.y;
    
    const Array1D &X = grid.X();
    for(unit_t i = lower.x; i <= upper.x; ++i)
    {
        new (jvert+i) Junction::List( Junction::Vert, X[i] );
    }
    
    const Array1D &Y = grid.Y();
    for(unit_t j=lower.y; j<= upper.y; ++j)
    {
        new (jhorz+j) Junction::List( Junction::Horz, Y[j] );
    }
    
}


Junctions:: ~Junctions() throw()
{
    for(unit_t i = upper.x; i >= lower.x; --i)
    {
        destruct(jvert+i);
    }
    
    for(unit_t j= upper.y; j >= lower.y; --j)
    {
        destruct(jhorz+j);
    }
    
    memory::kind<memory::global>::release_as(jlists, jcount);
}
