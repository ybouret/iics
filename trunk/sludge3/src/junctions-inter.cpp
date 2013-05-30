#include "junctions.hpp"


void Junctions:: inter( Bubble &bubble)
{
    assert(bubble.size>=3);
    bubble.flags = SLUDGE_INSIDE;
    
    //==========================================================================
    //
    // First pass: detect locations of tracers
    //
    //==========================================================================
    const Tracer *tr = bubble.root;
    for(size_t k=bubble.size;k>0;--k,tr=tr->next)
    {
        tr->flags = __Grid::Locate(grid, tr->pos, tr->coord);
        bubble.flags |= tr->flags;
    }
    
    //==========================================================================
    //
    // Second pass: create and count junctions
    //
    //==========================================================================
    assert(bubble.root == tr);
    for(size_t k=bubble.size;k>0;--k,tr=tr->next)
    {
        __intersect(bubble,tr);
    }
    
    
}


#include "yocto/code/utils.hpp"

void Junctions:: __intersect(const Bubble &bubble, const Tracer *u)
{
    assert(u!=NULL);
    assert(bubble.owns(u));
    assert(u->next);
    
    //==========================================================================
    //
    // copy coordinates : 81 possibilities
    //
    //==========================================================================
    const Vertex Ru = u->pos;
    const Coord  Cu = u->coord;
    const size_t Su = u->flags;
    
    const Tracer *v  = u->next;
    const Vertex  Rv = v->pos;
    const Coord   Cv = v->coord;
    const size_t  Sv = v->flags;
    
    //==========================================================================
    //
    // order by status: 45 possibilities left
    //
    //==========================================================================
    if(Sv<Su)
    {
        cswap_const(Ru,Rv);
        cswap_const(Cu,Cv);
        cswap_const(Su,Sv);
    }
    
    //==========================================================================
    //
    // study possibilities
    //
    //==========================================================================
    
    //--------------------------------------------------------------------------
    //
    // is it possibile to have an Horizontal junction ?
    //
    //--------------------------------------------------------------------------
    if( 0 == ( Su & SLUDGE_TOP_OR_BOTTOM) )
    {
        assert( Cu.y != SLUDGE_INVALID_COORD);
        //----------------------------------------------------------------------
        // possible if:
        // - v is at TOP or BOTTOM => automatic
        // - otherwise, y coordinates (must be valid) are different
        //----------------------------------------------------------------------
        if((0 != ( Sv & SLUDGE_TOP_OR_BOTTOM) ) ||
           (Cu.y != Cv.y) )
        {
            
#ifndef NDEBUG
            if( 0 ==(Sv&SLUDGE_TOP_OR_BOTTOM) )
            {
                assert(Cv.y!=SLUDGE_INVALID_COORD);
            }
#endif
            
            __interHorz(bubble, Ru, Cu, Rv);
        }
    }
    
    //--------------------------------------------------------------------------
    //
    // is it possibile to have an Vertical junction ?
    //
    //--------------------------------------------------------------------------
    if( 0 == ( Su & SLUDGE_LEFT_OR_RIGHT) )
    {
        assert(Cu.x != SLUDGE_INVALID_COORD);
        //----------------------------------------------------------------------
        // possible if:
        // - v is at LEFT or RIGHT
        // -- otherwise, x coordinates (must be valid) are differnet
        //----------------------------------------------------------------------
        if((0 != ( Sv & SLUDGE_LEFT_OR_RIGHT)) ||
           (Cu.x != Cv.x) )
        {
#ifndef NDEBUG
            if( 0 == (Sv&SLUDGE_LEFT_OR_RIGHT) )
            {
                assert(Cv.x!=SLUDGE_INVALID_COORD);
            }
#endif
            __interVert(bubble, Ru, Cu, Rv);
        }
    }
    
}

void Junctions:: __interHorz(const Bubble &bubble,
                             const Vertex &p,
                             const Coord  &P,
                             const Vertex &q)
{
    assert(P.y != SLUDGE_INVALID_COORD);
    const unit_t    j  = q.y>p.y ? P.y+1 : P.y;
    Junction::List &J  = Horz(j);
    const Real      y0 = J.level;
    const Real      x0 = p.x + (y0-p.y)*(q.x -p.x) /(q.y-p.y);
    J.append(x0);
}

void Junctions:: __interVert(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q)
{
    assert(P.x != SLUDGE_INVALID_COORD);
    const unit_t    i  = q.x>p.x ? P.x+1 : P.x;
    Junction::List &J  = Vert(i);
    const Real      x0 = J.level;
    const Real      y0 = p.y + (x0-p.x) * (q.y - p.y) / (q.x - p.x);
    J.append(y0);
}

