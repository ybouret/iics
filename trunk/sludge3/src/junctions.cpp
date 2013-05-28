#include "junctions.hpp"

////////////////////////////////////////////////////////////////////////////////
//
// JUNCTION
//
////////////////////////////////////////////////////////////////////////////////

Junction:: Junction( List &r, Real a) throw() :
prev(0),
next(0),
root(r),
value(a)
{
}

Junction:: ~Junction() throw()
{
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

void Junction::List::append( Real a)
{
    Junction *J = new Junction(*this,a);
    push_back(J);
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

void Junctions:: clear() throw()
{
    for(unit_t i = grid.upper.x; i >= grid.lower.x; --i)
        jvert[i].auto_delete();
    for(unit_t j= grid.upper.y; j >= grid.lower.y; --j)
        jhorz[j].auto_delete();
}


////////////////////////////////////////////////////////////////////////////////
//
// Junction Load algorithm
//
////////////////////////////////////////////////////////////////////////////////
bool Junctions:: load( const Bubble &bubble )
{
    assert(bubble.size>=3);
    const Tracer *tr = bubble.root;
    unsigned count = 0;
    for(size_t k=bubble.size;k>0;--k,tr=tr->next)
    {
        const Vertex &p = tr->pos;
        const Vertex &q = tr->next->pos;
        count += __load(bubble,p,q);
    }
    
    
    return count>0;
}

unsigned Junctions:: __load( const Bubble &bubble, const Vertex &p, const Vertex &q )
{
    Coord P;
    const int  ploc = __Grid::Locate(grid, p, P);
    const bool p_in = ploc == SLUDGE_INSIDE;
    
    Coord Q;
    const int  qloc = __Grid::Locate(grid, q, Q);
    const bool q_in = qloc == SLUDGE_INSIDE;
    
    if( p_in )
    {
        //----------------------------------------------------------------------
        // p is INSIDE
        //----------------------------------------------------------------------
        if(q_in)
        {
            //------------------------------------------------------------------
            //-- p and q are INSIDE
            //------------------------------------------------------------------
            return __loadJ(bubble, p, P, q, Q);
        }
        else
        {
            //------------------------------------------------------------------
            //-- p is INSIDE, q is OUTSIDE
            //------------------------------------------------------------------
            return __loadJ(bubble, p, P, q, Q);
        }
    }
    else
    {
        //----------------------------------------------------------------------
        // p is OUTSIDE
        //----------------------------------------------------------------------
        if(q_in)
        {
            //------------------------------------------------------------------
            // p is OUTSIDE, q is INSIDE
            //------------------------------------------------------------------
            return __loadJ(bubble, q, Q, p, P);
        }
        else
        {
            //------------------------------------------------------------------
            // p and q are OUTSIDE
            //------------------------------------------------------------------
            return 0;
        }
        
    }
}

unsigned Junctions:: __loadJ(const Bubble &bubble,
                             const Vertex &p,
                             const Coord  &P,
                             const Vertex &q,
                             const Coord  &Q)
{
    if( P.x == Q.x )
    {
        //======================================================================
        // SAME COLUMNS
        //======================================================================
        if( P.y == Q.y )
        {
            //------------------------------------------------------------------
            // SAME COLUMNS, SAME LINES => nothing
            //------------------------------------------------------------------
            //std::cerr << "same case" << std::endl;
            return 0;
        }
        else
        {
            //------------------------------------------------------------------
            // SAME COLUMNS, DIFFERENT LINES => will cut an horizontal axis
            //------------------------------------------------------------------
            //assert(Q.y == P.y-1 || Q.y == P.y+1);
            __loadHorz(bubble, p,P,q);
            return 1;
        }
    }
    else
    {
        //assert(Q.x==P.x-1||Q.x==P.x+1);
        //======================================================================
        // DIFFERENT COLUMNS
        //======================================================================
        if( P.y == Q.y )
        {
            //------------------------------------------------------------------
            // DIFFERENT COLUMNS, SAME LINES => will cut a vertical axis
            //------------------------------------------------------------------
            __loadVert(bubble, p, P, q);
            return 1;
        }
        else
        {
            //------------------------------------------------------------------
            // DIFFERENT COLUMNS, DIFFERENT LINES => will cut two axis
            //------------------------------------------------------------------
            __loadHorz(bubble, p, P, q);
            __loadVert(bubble, p, P, q);
            return 2;
        }
    }
    
}


void Junctions:: __loadHorz(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q)
{
    const unit_t    j  = q.y>p.y ? P.y+1 : P.y;
    Junction::List &J  = Horz(j);
    const Real      y0 = J.level;
    const Real      x0 = p.x + (y0-p.y)*(q.x -p.x) /(q.y-p.y);
    J.append(x0);
}

void Junctions:: __loadVert(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q)
{
    const unit_t    i  = q.x>p.x ? P.x+1 : P.x;
    Junction::List &J  = Vert(i);
    const Real      x0 = J.level;
    const Real      y0 = p.y + (x0-p.x) * (q.y - p.y) / (q.x - p.x);
    J.append(y0);
}

////////////////////////////////////////////////////////////////////////////////
//
// Junctions sorting
//
////////////////////////////////////////////////////////////////////////////////
#include "yocto/core/merge-sort.hpp"
#include "yocto/comparator.hpp"

static inline
int __compareJ( const Junction *lhs, const Junction *rhs, void *)
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


