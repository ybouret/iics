#include "junctions.hpp"


size_t Junctions:: inter( const Bubble &bubble)
{
    assert(bubble.size>=3);
    size_t ans = SLUDGE_INSIDE;
    
    //==========================================================================
    //
    // First pass: detect locations of tracers
    //
    //==========================================================================
    const Tracer *tr = bubble.root;
    for(size_t k=bubble.size;k>0;--k,tr=tr->next)
    {
        tr->flags = __Grid::Locate(grid, tr->pos, tr->coord);
        ans |= tr->flags;
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
    
    
    return ans;
}



void Junctions:: __intersect(const Bubble &bubble, const Tracer *u)
{
    assert(u!=NULL);
    assert(bubble.owns(u));
    assert(u->next);
    
    //==========================================================================
    //
    // copy coordinates
    //
    //==========================================================================
    const Vertex Ru = u->pos;
    const Coord  Cu = u->coord;
    const size_t Su = u->flags;
    
    const Tracer *v = u->next;
    const Vertex  Rv = v->pos;
    const Coord   Cv = v->coord;
    const size_t  Sv = v->flags;
    
    
    
    
}