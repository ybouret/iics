#include "../segmenter.hpp"
#include "yocto/core/merge-sort.hpp"

static inline
int __compare_horz( const Junction *lhs, const Junction *rhs, void * ) throw()
{
    return __compare<Real>(lhs->vertex.x,rhs->vertex.x);
}

static inline
int __compare_vert( const Junction *lhs, const Junction *rhs, void * ) throw()
{
    return __compare<Real>(lhs->vertex.y,rhs->vertex.y);
}

void Segmenter:: SortHorz()
{
    for( unit_t j=Y.lower;j<=Y.upper;++j)
    {
        core::list_of<Junction> &js = Horz(j);
        core::merging<Junction>::sort( js, __compare_horz,0);
    }
}

void Segmenter:: SortVert()
{
    for(unit_t i=X.lower;i<=X.upper;++i)
    {
		core::list_of<Junction> &js = Vert(i);
        core::merging<Junction>::sort( js, __compare_vert,0);
    }
}

void Segmenter:: remove_vertical_junctions_below( const Real ylim )
{
    
    for(unit_t i=X.lower;i<=X.upper;++i)
    {
        Junctions &jvert = Vert(i);
        core::list_of<Junction> tmp;
        while( jvert.size )
        {
            Junction *J = jvert.pop_front();
            if( J->vertex.y < ylim )
            {
                duplicates.push_back(J);
            }
            else
            {
                tmp.push_back(J);
            }
        }
        jvert.swap_with(tmp);
    }
}

void Segmenter:: process( const Bubbles &bubbles )
{
    assert(hseg!=NULL);
    assert(vseg!=NULL);
    
    //--------------------------------------------------------------------------
    // reset
    //--------------------------------------------------------------------------
    for( size_t i=segcount;i>0;--i)
        segments[i]->empty();
    markers.empty();
    duplicates.empty();
    
    //--------------------------------------------------------------------------
    // find all junctions
    //--------------------------------------------------------------------------
    for( const Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next)
    {
        process_bubble(bubble, bubbles.pbc.up);
    }
    
    //--------------------------------------------------------------------------
    // sort junctions
    //--------------------------------------------------------------------------
    SortHorz();
    SortVert();
    
    //--------------------------------------------------------------------------
    // locate Horizontal junctions (for building B field + gradient eval)
    //--------------------------------------------------------------------------
    for( unit_t j=Y.lower;j<=Y.upper;++j)
    {
        for( Junction *J = Horz(j).head; J; J=J->next)
        {
            locate_value( J->vertex.x, X, J->klo, J->khi);
        }
    }
    
    //--------------------------------------------------------------------------
    // locate Vertical junctions (for gradient eval)
    //--------------------------------------------------------------------------
    for( unit_t i=X.lower;i<=X.upper;++i)
    {
        for( Junction *J = Vert(i).head; J; J=J->next)
        {
            locate_value( J->vertex.y, Y, J->klo, J->khi);
        }
    }
    
    
}

void Segmenter:: process_bubble( const Bubble *bubble , const Real half)
{
    assert(bubble);
    
    //--------------------------------------------------------------------------
    // outer loop on bubble spots
    //--------------------------------------------------------------------------
    for( const Spot *spot=bubble->spots.head;spot;spot=spot->next)
    {
        process_spot(spot,half);
    }
    
}

void Segmenter:: process_spot( const Spot *spot, const Real half )
{
    const Tracer *tracer = spot->handle;
    const Vertex  vertex = tracer->vertex;
    assert( vertex.x >= X[X.lower]);
    assert( vertex.x <  X[X.upper]);
    assert( vertex.y >= Y[Y.lower]);
    assert( vertex.y <  Y[Y.upper]);
    
    //--------------------------------------------------------------------------
    // locate the tracer
    //--------------------------------------------------------------------------
    locate_vertex(vertex, spot->klo, spot->kup);
    
    //--------------------------------------------------------------------------
    // compute the absolute next coordinate
    //--------------------------------------------------------------------------
    Vertex target = vertex + tracer->edge;
    
    //--------------------------------------------------------------------------
    // find junctions
    //--------------------------------------------------------------------------
    compute_junctions(spot, vertex, target, tracer->next);
    const Tracer *prec = tracer->prev;
    if( ! prec->is_spot || Fabs(prec->vertex.y - vertex.y) >= half )
    {
        target = vertex - prec->edge;
        compute_junctions(spot, vertex, target, prec);
    }
}

typedef int OutCode;

const OutCode INSIDE = 0; // 0000
const OutCode LEFT   = 1; // 0001
const OutCode RIGHT  = 2; // 0010
const OutCode BOTTOM = 4; // 0100
const OutCode TOP    = 8; // 1000


// Compute the bit code for a point (x, y) using the clip rectangle
// bounded diagonally by (xmin, ymin), and (xmax, ymax)

// ASSUME THAT xmax, xmin, ymax and ymin are global constants.
static inline
OutCode ComputeOutCode(const Real x, const Real y, const Real xmin, const Real xmax, const Real ymin, const Real ymax) throw()
{
    OutCode code;
    
    code = INSIDE;           // initialised as being inside of clip window
    
    if (x < xmin)            // to the left of clip window
        code |= LEFT;
    else if (x >= xmax)      // to the right of clip window
        code |= RIGHT;
    if (y < ymin)            // below the clip window
        code |= BOTTOM;
    else if (y >= ymax)      // above the clip window
        code |= TOP;
    
    return code;
}

static inline Real dist_sq( const Vertex &lhs, const Vertex &rhs ) throw()
{
    const Vertex d(lhs,rhs);
    return d.norm2();
}

static inline
void update_jnext( Tracer *p, const Junction *J )
{
    if(
       ( p->jnext == 0 ) ||
       ( dist_sq(p->vertex,J->vertex) < dist_sq(p->vertex, p->jnext->vertex) )
       )
    {
        p->jnext = J;
    }
    
}

static inline
void update_jprev( Tracer *p, const Junction *J )
{
    if(
       ( p->jprev == 0 ) ||
       ( dist_sq(p->vertex,J->vertex) < dist_sq(p->vertex, p->jprev->vertex) )
       )
    {
        p->jprev = J;
    }
    
}

static inline void FinalizeJunction( Junction *J, const Tracer *source, const Tracer *target ) throw()
{
    assert( target == source->next || target == source->prev );
    const Real  s_weight = (1-J->alpha);
    const Real  t_weight = J->alpha;
    J->bubble            = source->bubble;
    J->curvature         = s_weight*source->curvature + t_weight * target->curvature;
    const Real s_angle   = source->angle;
    const Real t_angle   = target->angle;
    const Real j_angle   = s_weight * s_angle + t_weight * t_angle;
    J->n.x = Cos(j_angle);
    J->n.y = Sin(j_angle);
    
    if( target == source->next )
    {
        update_jnext( (Tracer *)source, J);
        update_jprev( (Tracer *)target, J);
        return;
    }
    
    if( target == source->prev )
    {
        update_jprev( (Tracer *)source, J);
        update_jnext( (Tracer *)target, J);
    }
}

static inline
OutCode ComputeOutCode( const Vertex &r, const Vertex &lo, const Vertex &up) throw()
{
    return ComputeOutCode(r.x, r.y, lo.x, up.x, lo.y, up.y);
}


void Segmenter:: compute_junctions(const Spot   *spot,
                                   const Vertex &self,
                                   const Vertex &other,
                                   const Tracer *to)
{
    const Coord   klo = spot->klo;
    const Coord   kup = spot->kup;
    
    const Vertex lo( X[klo.x], Y[klo.y]);
    const Vertex up( X[kup.x], Y[kup.y]);
    assert( INSIDE == ComputeOutCode(self,lo,up));
    
    const OutCode other_code = ComputeOutCode(other, lo, up);
    
    if( other_code & LEFT)
    {
        assert( other.x < self.x);
        Segment  &seg = Vert(klo.x);
        Junction *J   = seg.append();
        J->vertex.x   = seg.value;
        J->alpha      = ( J->vertex.x - self.x)/(other.x - self.x);
        J->vertex.y   = self.y + (J->alpha)*(other.y - self.y);
        FinalizeJunction(J,spot->handle,to);
    }
    
    if( other_code & RIGHT)
    {
        assert( other.x > self.x);
        Segment  &seg = Vert(kup.x);
        Junction *J   = seg.append();
        J->vertex.x   = seg.value;
        J->alpha      = ( J->vertex.x - self.x)/(other.x - self.x);
        J->vertex.y   = self.y + (J->alpha)*(other.y - self.y);
        FinalizeJunction(J,spot->handle,to);
    }
    
    if( other_code & BOTTOM )
    {
        assert( other.y < self.y );
        Segment  &seg = Horz(klo.y);
        Junction *J   = seg.append();
        J->vertex.y   = seg.value;
        J->alpha      = (J->vertex.y - self.y)/(other.y-self.y);
        J->vertex.x   = self.x + (J->alpha)*(other.x-self.x);
        FinalizeJunction(J,spot->handle,to);
    }
    
    
    if( other_code & TOP )
    {
        assert( other.y > self.y );
        Segment  &seg = Horz(kup.y);
        Junction *J   = seg.append();
        J->vertex.y   = seg.value;
        J->alpha      = (J->vertex.y - self.y)/(other.y-self.y);
        J->vertex.x   = self.x + (J->alpha)*(other.x-self.x);
        FinalizeJunction(J,spot->handle,to);
    }
    
    
}
