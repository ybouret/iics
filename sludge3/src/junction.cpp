#include "junction.hpp"

////////////////////////////////////////////////////////////////////////////////
//
// JUNCTION
//
////////////////////////////////////////////////////////////////////////////////

Junction:: Junction( List &r, Real a, const Bubble *b) throw() :
prev(0),
next(0),
root(r),
value(a),
owner(b),
C(0),
inside(false),
lower( SLUDGE_INVALID_COORD ),
upper( SLUDGE_INVALID_COORD ),
b_pos( Bubble::IsInvalid )
{
    assert(owner);
}

Junction:: ~Junction() throw()
{
}


Vertex Junction:: get(void) const throw()
{
    switch(root.type)
    {
        case Horz:
            return  Vertex(value,root.level);
            
        case Vert:
            return Vertex(root.level,value);
    }
}


////////////////////////////////////////////////////////////////////////////////
//
// JUNCTION LIST
//
////////////////////////////////////////////////////////////////////////////////
Junction::List:: List( Type t, const Real &v) throw() :
type(t),
level(v)
{
}

Junction::List:: ~List() throw() { auto_delete(); }

Junction *Junction::List::append( Real value, const Bubble *owner)
{
    assert(owner);
    Junction *J = new Junction(*this,value,owner);
    push_back(J);
    return J;
}

const Junction * Junction:: List::after(const Real value) const throw()
{
    for(const Junction *J=head;J;J=J->next)
    {
        if(J->value>=value)
            return J;
    }
    return 0;
}


const Junction * Junction:: List::before(const Real value) const throw()
{
    for(const Junction *J=tail;J;J=J->prev)
    {
        if(J->value<=value)
            return J;
    }
    return 0;
}



