#include "segmenter.hpp"

Segmenter:: Segmenter( const Grid &grid ) :
X(  grid.X()  ),
Y(  grid.Y()  ),
dX( grid.dX() ),
dY( grid.dY() ),
j_cache(),
junctions(j_cache),
s_cache(),
horizontal(0),
vertical(0),
segments()
{
}

Segmenter::~Segmenter() throw()
{
    
}

void Segmenter:: clear() throw()
{
    junctions.empty();
    for( size_t i=segments.size(); i>0; --i ) segments[i].empty();
    
}

void Segmenter:: allocate_segments()
{
    const size_t nx = X.width;
    const size_t ny = Y.width;
    const size_t n  = nx+ny;
    
    horizontal = vertical = 0;
    clear();
    {
        const Segment::List seg(s_cache);
        segments.free();
        segments.make(n,seg);
    }
    horizontal = &segments[1]    - Y.lower;
    vertical   = &segments[1+ny] - X.lower;
}

//==============================================================================
// bissection to locate where a point is
//==============================================================================
static unit_t __locate_point(Real           x, 
                             const Array1D &a)
{
    unit_t ilo = a.lower;
    if( x <= a[ilo] )
        return ilo;
    
    unit_t ihi = a.upper;
    assert(ihi-ilo>=1);
    if( x >= a[ihi] )
        return ihi -1;
    
    while(ihi-ilo>1)
    {
        const int mid = (ilo+ihi)>>1;
        if( x < a[mid] )
        {
            ihi = mid;
        }
        else 
        {
            ilo = mid;
        }
    }
    assert(a[ilo]<=x);
    assert(x<a[ihi]);
    return ilo;
}


void Segmenter:: process_tracer(Tracer *p)
{
    assert(segments.size()>0);
    assert(p!=NULL);
    const Vertex P = p->vertex;
    Bubble      *bubble = p->bubble; assert(bubble);
    assert(P.y >= Y[Y.lower] );
    assert(P.y <  Y[Y.upper] );
    assert(P.x >= X[X.lower] );
    assert(P.x <= X[X.upper] );
    
    //--------------------------------------------------------------------------
    // find position on grid
    //--------------------------------------------------------------------------
    const unit_t i = p->gpos.x = __locate_point( P.x, X );
    const unit_t j = p->gpos.y = __locate_point( P.y, Y );
    
    //--------------------------------------------------------------------------
    // find grid boundaries
    //--------------------------------------------------------------------------
    const Vertex vmin( X[i],   Y[j]   );
    const Vertex vmax( X[i+1], Y[j+1] );
    
    //--------------------------------------------------------------------------
    // compute bilinear interpolation coefficients
    //--------------------------------------------------------------------------
    p->bw.x = (P.x-vmin.x)/dX[i];
    p->bw.y = (P.y-vmin.y)/dY[j];
    
    //--------------------------------------------------------------------------
    // find potential intersections
    //--------------------------------------------------------------------------
    {
        const Vertex Q = P + p->edge;
        find_junctions(P, Q, vmin, vmax, p);
    }

    if( !p->prev->is_spot )
    {
        const Vertex Q = P - (p->prev->edge);
        find_junctions(P, Q, vmin, vmax, p);
    }
    
}


#define __NEW_INTER(SEGMENT,INDEX,COORD,ARR)  do  { \
Junction *I = junctions.append(); \
I->vertex.x = Ix;   \
I->vertex.y = Iy;   \
I->bubble   = bubble; \
SEGMENT[INDEX].append()->handle = I; \
I->up = 1 + (I->lo = __locate_point( I->vertex.COORD, ARR )); \
} while(false)



void Segmenter:: find_junctions(const Vertex &P, 
                                const Vertex &Q, 
                                const Vertex &vmin, 
                                const Vertex &vmax, 
                                Tracer       *p )
{
    assert(p!=NULL);
    const unit_t i = p->gpos.x;
    const unit_t j = p->gpos.y;
    const unit_t i1 = i+1;
    const unit_t j1 = j+1;
    Bubble *bubble  = p->bubble; assert(bubble);
    //--------------------------------------------------------------------------
    // simplified Cohen-Sutherland
    //--------------------------------------------------------------------------
    if( Q.x < vmin.x )
    {
        assert( Q.x < P.x );
        const Real Iy = P.y + ( vmin.x - P.x )  * (Q.y - P.y ) / (Q.x - P.x );
        const Real Ix = vmin.x;
        __NEW_INTER(vertical,i,y,Y);
    }
    
    // special right case
    if( P.x < vmax.x && Q.x >= vmax.x )
    {
        assert(Q.x>P.x);
        const Real Iy = P.y + ( vmax.x - P.x )  * (Q.y - P.y ) / (Q.x - P.x );
        const Real Ix = vmax.x;
        __NEW_INTER(vertical,i1,y,Y);
    }
    
    if( Q.y < vmin.y )
    {
        assert( Q.y < P.y );
        const Real Ix = P.x + ( vmin.y - P.y)* (Q.x - P.x) /(Q.y - P.y );
        const Real Iy = vmin.y;
        __NEW_INTER(horizontal,j,x,X);
    }
    
    if( Q.y >= vmax.y )
    {
        assert( Q.y > P.y );
        const Real Ix = P.x + ( vmax.y - P.y)* (Q.x - P.x) / (Q.y - P.y ) ;
        const Real Iy = vmax.y;
        __NEW_INTER(horizontal,j1,x,X);
    }
    
}



void Segmenter:: process_bubble(Bubble *bubble)
{
    assert(bubble);
    bubble->markers.empty();
    for( Spot *spot = bubble->spots.head; spot; spot=spot->next )
    {
        process_tracer( spot->handle );
    }
}

void Segmenter:: process_bubbles( Bubbles &bubbles )
{
    clear();
    for( Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next)
    {
        process_bubble(bubble);
    }
    sort_segments();
}

#include "yocto/core/merge-sort.hpp"

static inline int __compare_segments( const Segment *lhs, const Segment *rhs, void *) throw()
{
    return lhs->handle->lo - rhs->handle->lo;
}

void Segmenter:: sort_segments()
{
    for( size_t i=segments.size(); i>0; --i ) 
    {
        core::merging<Segment>::sort( segments[i], __compare_segments, 0 );
    }
}

#include "yocto/ios/ocstream.hpp"
void Segmenter:: save_junctions( const string &filename ) const
{
    ios::ocstream fp( filename, false );
    for( const Junction *J = junctions.head; J; J=J->next )
    {
        fp("%g %g\n", J->vertex.x, J->vertex.y );
    }
}
