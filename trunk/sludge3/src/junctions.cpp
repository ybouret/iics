#include "junctions.hpp"

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
upper( SLUDGE_INVALID_COORD )
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

////////////////////////////////////////////////////////////////////////////////
//
// JUNCTIONS
//
////////////////////////////////////////////////////////////////////////////////
#include "yocto/memory/global.hpp"

Junctions:: Junctions(Grid &g) :
grid(g),
num_lists(grid.width.x + grid.width.y),
jcount(num_lists),
jlists( memory::kind<memory::global>::acquire_as<Junction::List>(jcount) ),
jvert(jlists),
jhorz(jvert+grid.width.x)
{
    jvert -= grid.lower.x;
    jhorz -= grid.lower.y;
    
    const Array1D &X = grid.X();
    for(unit_t i = grid.lower.x; i <= grid.upper.x; ++i)
    {
        new (jvert+i) Junction::List( Junction::Vert, X[i] );
    }
    
    const Array1D &Y = grid.Y();
    for(unit_t j=grid.lower.y; j<= grid.upper.y; ++j)
    {
        new (jhorz+j) Junction::List( Junction::Horz, Y[j] );
    }
    
}


Junctions:: ~Junctions() throw()
{
    for(unit_t i = grid.upper.x; i >= grid.lower.x; --i)
    {
        destruct(jvert+i);
    }
    
    for(unit_t j= grid.upper.y; j >= grid.lower.y; --j)
    {
        destruct(jhorz+j);
    }
    
    memory::kind<memory::global>::release_as(jlists, jcount);
}

Junction::List & Junctions:: Vert( unit_t i ) throw()
{
    assert(i>=grid.lower.x);
    assert(i<=grid.upper.x);
    return jvert[i];
}

Junction::List & Junctions:: Horz( unit_t j ) throw()
{
    assert(j>=grid.lower.y);
    assert(j<=grid.upper.y);
    return jhorz[j];
}

const Junction::List & Junctions:: Vert( unit_t i ) const throw()
{
    assert(i>=grid.lower.x);
    assert(i<=grid.upper.x);
    return jvert[i];
}

const Junction::List & Junctions:: Horz( unit_t j ) const throw()
{
    assert(j>=grid.lower.y);
    assert(j<=grid.upper.y);
    return jhorz[j];
}


void Junctions:: clear() throw()
{
    for(unit_t i = grid.upper.x; i >= grid.lower.x; --i)
        jvert[i].auto_delete();
    for(unit_t j= grid.upper.y; j >= grid.lower.y; --j)
        jhorz[j].auto_delete();
}



////////////////////////////////////////////////////////////////////////////////
//
// Junctions sorting
//
////////////////////////////////////////////////////////////////////////////////
#include "yocto/core/merge-sort.hpp"
#include "yocto/comparator.hpp"

static inline
int __compareJ( const Junction *lhs, const Junction *rhs, void *) throw()
{
    assert(lhs);
    assert(rhs);
    return __compare<Real>(lhs->value,rhs->value);
}

void Junctions:: sort()
{
    for( size_t i=0; i<num_lists; ++i)
    {
        core::merging<Junction>::sort<core::list_of>(jlists[i],__compareJ,(void*)0);
    }
}

////////////////////////////////////////////////////////////////////////////////
//
// Junctions I/O
//
////////////////////////////////////////////////////////////////////////////////

#include "yocto/ios/ocstream.hpp"

void Junctions:: save_dat( const string &fn) const
{
    ios::ocstream fp( fn, false);
    
    for(unit_t i = grid.upper.x; i >= grid.lower.x; --i)
    {
        const Junction::List &JL = jvert[i];
        for( const Junction *J = JL.head; J; J=J->next)
        {
            fp("%g %g\n", J->root.level, J->value);
        }
    }
    
    for(unit_t j= grid.upper.y; j >= grid.lower.y; --j)
    {
        const Junction::List &JL = jhorz[j];
        for( const Junction *J = JL.head; J; J=J->next)
        {
            fp("%g %g\n", J->value, J->root.level);
        }
        
    }
}

void Junctions:: save_t( const string &fn ) const
{
    ios::ocstream fp( fn, false);
    for(size_t i=0; i < num_lists; ++i )
    {
        const Junction::List &JL = jlists[i];
        for( const Junction *J = JL.head;J;J=J->next)
        {
            const Vertex p = J->get();
            fp("%g %g\n", p.x, p.y);
            const Vertex q = p + J->owner->lambda * J->t;
            fp("%g %g\n", q.x, q.y);
            fp("\n");
        }
    }
}

void Junctions:: save_n( const string &fn ) const
{
    ios::ocstream fp( fn, false);
    for(size_t i=0; i < num_lists; ++i )
    {
        const Junction::List &JL = jlists[i];
        for( const Junction *J = JL.head;J;J=J->next)
        {
            const Vertex p = J->get();
            fp("%g %g\n", p.x, p.y);
            const Vertex q = p + (J->owner->lambda * J->C ) * J->n;
            fp("%g %g\n", q.x, q.y);
            fp("\n");
        }
    }
}


