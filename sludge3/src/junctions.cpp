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

Junction * Junction::List::append( Real a)
{
    Junction *J = new Junction(*this,a);
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
jcount( grid.width.x + grid.width.y ),
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
            return __loadBothInside(bubble, p, P, q, Q);
        }
        else
        {
            //------------------------------------------------------------------
            //-- p is INSIDE, q is OUTSIDE
            //------------------------------------------------------------------
            return __loadOneOutside(bubble, p, P, q);
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
            return __loadOneOutside(bubble, q, Q, p);
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

unsigned Junctions:: __loadBothInside(const Bubble &bubble,
                                      const Vertex &p,
                                      const Coord  &P,
                                      const Vertex &q,
                                      const Coord  &Q)
{
    return 0;
}

unsigned Junctions:: __loadOneOutside(const Bubble &bubble,
                                      const Vertex &p,
                                      const Coord  &P,
                                      const Vertex &q)
{
    return 0;
}
