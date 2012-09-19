#include "segmenter.hpp"

void Segmenter:: process( const Bubbles &bubbles )
{
    assert(hseg!=NULL);
    assert(vseg!=NULL);
    for( size_t i=segcount;i>0;--i) segments[i]->empty();
    for( const Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next)
    {
        process_bubble(bubble);
    }
}

void Segmenter:: process_bubble( const Bubble *bubble )
{
    assert(bubble);
    
    //--------------------------------------------------------------------------
    // outer loop on bubble spots
    //--------------------------------------------------------------------------
    for( const Spot *spot=bubble->spots.head;spot;spot=spot->next)
    {
        process_spot(spot);
    }
    
}

void Segmenter:: process_spot( const Spot *spot )
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
    const Vertex target = vertex + tracer->edge;
    
    //--------------------------------------------------------------------------
    // find junctions
    //--------------------------------------------------------------------------
    compute_junctions(spot, vertex, target);
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

static inline
OutCode ComputeOutCode( const Vertex &r, const Vertex &lo, const Vertex &up)
{
    return ComputeOutCode(r.x, r.y, lo.x, up.x, lo.y, up.y);
}

void Segmenter:: compute_junctions(const Spot   *spot,
                                   const Vertex &self,
                                   const Vertex &other )
{
    const Bubble *bubble = spot->handle->bubble;
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
        J->pos.x      = seg.value;
        J->pos.y      = self.y + ( J->pos.x - self.x)*(other.y - self.y)/(other.x - self.x) ;
        J->bubble     = bubble;
    }
    
    if( other_code & RIGHT)
    {
        assert( other.x > self.x);
        Segment  &seg = Vert(kup.x);
        Junction *J   = seg.append();
        J->pos.x      = seg.value;
        J->pos.y      = self.y + ( J->pos.x - self.x)*(other.y - self.y)/(other.x - self.x);
        J->bubble     = bubble;
    }

    if( other_code & BOTTOM )
    {
        assert( other.y < self.y );
        Segment  &seg = Horz(klo.y);
        Junction *J   = seg.append();
        J->pos.y      = seg.value;
        J->pos.x      = self.x + (J->pos.y - self.y)*(other.x-self.x)/(other.y-self.y);
        J->bubble     = bubble;
    }

    
    if( other_code & TOP )
    {
        assert( other.y > self.y );
        Segment  &seg = Horz(kup.y);
        Junction *J   = seg.append();
        J->pos.y      = seg.value;
        J->pos.x      = self.x + (J->pos.y - self.y)*(other.x-self.x)/(other.y-self.y);
        J->bubble     = bubble;
    }
    
      
    
    
}
