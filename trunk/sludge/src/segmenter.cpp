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
    
}

void Segmenter:: horizontal_pbc( unit_t y_lower, unit_t y_upper )
{
    assert(y_lower<y_upper);
    assert(y_lower>=Y.lower);
    assert(y_upper<=Y.upper);
    Segment::List &S_lo = horizontal[y_lower];
    Segment::List &S_up = horizontal[y_upper];
    
    const Real y_lo = Y[y_lower];
    const Real y_up = Y[y_upper];
    
    Segment::List mirror_lo( s_cache );
    for( const Segment *s = S_lo.head; s; s=s->next )
    {
        const Junction *J = s->handle;
        Junction       *K = junctions.append();
        K->copy(J);
        K->vertex.y = y_up;
        mirror_lo.attach(K);
    }
    
    Segment::List mirror_up( s_cache );
    for( const Segment *s = S_up.head; s; s=s->next )
    {
        const Junction *J = s->handle;
        Junction       *K = junctions.append();
        K->copy(J);
        K->vertex.y = y_lo;
        mirror_up.attach(K);
    }
    
    S_up.merge_back(mirror_lo);
    S_lo.merge_back(mirror_up);
    
    
    
}


#include "yocto/core/merge-sort.hpp"
#include "yocto/comparator.hpp"
static inline int __compare_horz( const Segment *lhs, const Segment *rhs, void *) throw()
{
    return __compare( lhs->handle->vertex.x, rhs->handle->vertex.x );
}

static inline int __compare_vert( const Segment *lhs, const Segment *rhs, void *) throw()
{
    return __compare( lhs->handle->vertex.y, rhs->handle->vertex.y );
}

void Segmenter:: sort_segments()
{
    for( unit_t j=Y.lower;j<=Y.upper;++j)
        core::merging<Segment>::sort( horizontal[j], __compare_horz, 0);
    
    
    for( unit_t i=X.lower;i<=X.upper;++i)
        core::merging<Segment>::sort( vertical[i], __compare_vert, 0);
    
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
