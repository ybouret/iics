#include "junctions.hpp"


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
    /*
    for(unit_t i = grid.upper.x; i >= grid.lower.x; --i)
        jvert[i].auto_delete();
    for(unit_t j= grid.upper.y; j >= grid.lower.y; --j)
        jhorz[j].auto_delete();
     */
    for( size_t i=0; i<num_lists; ++i)
    {
        jlists[i].auto_delete();
    }
}


size_t Junctions:: count_all() const throw()
{
    size_t ans = 0;
    for( size_t i=0; i<num_lists; ++i)
    {
        ans += jlists[i].size;
    }
    return ans;
}

void  Junctions:: to_curve( array<Real> &cx, array<Real> &cy ) const throw()
{
    const size_t nj = count_all();
    assert(cx.size()==nj);
    assert(cy.size()==nj);
    size_t k = 0;
    
    for( size_t i=0; i<num_lists; ++i)
    {
        for( const Junction *J = jlists[i].head;J;J=J->next)
        {
            const Vertex v = J->get();
            ++k;
            cx[k] = v.x;
            cy[k] = v.y;
        }
    }
    assert(nj==k);
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

void Junctions:: save_all(const string &pfx) const
{
    {
        const string fn = pfx + ".dat";
        save_dat(fn);
    }
    
    {
        const string fn = pfx + "_t.dat";
        save_t(fn);
    }
    
    {
        const string fn = pfx  + "_n.dat";
        save_n(fn);
    }
    
}


