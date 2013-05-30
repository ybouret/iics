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
    
    const size_t Su = u->flags;
    
    const Tracer *v  = u->next;
    const size_t  Sv = v->flags;
    
    //==========================================================================
    //
    // order by status: 45 possibilities left
    //
    //==========================================================================
    if(Sv<Su)
    {
        cswap(u,v);
        cswap_const(Su,Sv);
    }
    assert(u->flags<=v->flags);
    
    const Vertex  Ru = u->pos;
    const Coord   Cu = u->coord;
    const Vertex  Rv = v->pos;
    const Coord   Cv = v->coord;
    
    //==========================================================================
    //
    // study possibilities: only a few left !
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
            Real      alpha = -1;
            Junction *J     = __interHorz(bubble, Ru, Cu, Rv, alpha);
            __updateJunction(J,alpha,u,v);
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
            Real     alpha =-1;
            Junction *J    = __interVert(bubble, Ru, Cu, Rv,alpha);
            __updateJunction(J,alpha,u,v);
        }
    }
    
}

Junction *Junctions:: __interHorz(const Bubble &bubble,
                                  const Vertex &p,
                                  const Coord  &P,
                                  const Vertex &q,
                                  Real         &alpha)
{
    assert(P.y != SLUDGE_INVALID_COORD);
    const unit_t    j     = q.y>p.y ? P.y+1 : P.y;
    Junction::List &JL    = Horz(j);
    const Real      y0    = JL.level;
    alpha                 = clamp<Real>(0,(y0-p.y)/(q.y - p.y ),1);
    const Real      x0    = p.x + (q.x -p.x)*alpha;
    Junction       *J     = JL.append(x0,&bubble);
    const Array1D  &X     = grid.X();
    if(J->value>=X[X.lower]&&J->value<=X[X.upper])
    {
        J->inside = true;
        J->lower  = __Grid::FindLower(X, J->value);
        J->upper  = J->lower+1;
    }
    return J;
}

Junction * Junctions:: __interVert(const Bubble &bubble,
                                   const Vertex &p,
                                   const Coord  &P,
                                   const Vertex &q,
                                   Real         &alpha)
{
    assert(P.x != SLUDGE_INVALID_COORD);
    const unit_t    i  = q.x>p.x ? P.x+1 : P.x;
    Junction::List &JL = Vert(i);
    const Real      x0 = JL.level;
    alpha              = clamp<Real>(0, (x0-p.x)/(q.x - p.x), 1);
    const Real      y0 = p.y +  alpha * (q.y - p.y);
    Junction       *J  = JL.append(y0,&bubble);
    const Array1D  &Y  = grid.Y();
    if(J->value>=Y[Y.lower] && J->value<=Y[Y.upper])
    {
        J->inside = true;
        J->lower  = __Grid::FindLower(Y,J->value);
        J->upper  = J->upper+1;
    }
    return J;
}

void Junctions:: __updateJunction( Junction *J, const Real alpha, const Tracer *u, const Tracer *v)
{
    assert(J);
    assert(u);
    assert(v);
    const Real U = (1-alpha); // weight for u
    const Real V = alpha;     // weight for v
    
    //--------------------------------------------------------------------------
    //-- average curvature
    //--------------------------------------------------------------------------
    J->C = U * u->C + V * v->C;
    J->pressure = J->owner->pressure - J->owner->gamma * J->C;
    
    //--------------------------------------------------------------------------
    //-- compute angle
    //--------------------------------------------------------------------------
    const Real theta = Vertex::angle_of(u->t, v->t);
    J->t = u->t.rotated_by( V * theta );
    J->t.normalize();
    J->n.x = - J->t.y;
    J->n.y =   J->t.x;
}


void Junctions:: load( Bubbles &bubbles )
{
    clear();
    for( Bubble *b = bubbles.head; b; b=b->next )
    {
        inter( *b );
    }
    sort();
}

