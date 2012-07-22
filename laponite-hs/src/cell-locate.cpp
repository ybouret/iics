#include "cell.hpp"

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

//==============================================================================
//
// Cohen-Sutherland algorithm
//
//==============================================================================


//==============================================================================
// encode location
//==============================================================================

static const int INSIDE   = 0;  // 0000
static const int LEFT     = 1;  // 0001
static const int RIGHT    = 2;  // 0010
static const int BOTTOM   = 4;  // 0100
static const int TOP      = 8;  // 1000
static const int TEST_ALL = LEFT | RIGHT | BOTTOM | TOP;


static inline 
void __OutMask( int &flag, int no ) throw()
{
    flag &= ~no;
}


static inline 
int __OutCode(const V2D &q,
              const V2D &vmin, 
              const V2D &vmax
              ) throw()
{
    int ans = INSIDE;
    
    if( (q.x<vmin.x) )
    {
        ans |= LEFT;
    }
    
    if( (q.x >= vmax.x) )
    {
        ans |= RIGHT;
    }
    
    if(  (q.y< vmin.y) )
    {
        ans |= BOTTOM;
    }
    
    if( q.y >= vmax.y ) 
    {
        ans |= TOP;
    }
    
    return ans;
}

//==============================================================================
// successive clipping
//==============================================================================

#define __NEW_INTER(SEGMENT,INDEX,COORD,ARR)  do  { \
Intersection *I = inter.append(); \
I->vertex.x = Ix;   \
I->vertex.y = Iy;   \
I->bubble = bubble; \
SEGMENT[INDEX].append()->inter = I; \
I->up = 1 + (I->lo = __locate_point( I->vertex.COORD, ARR )); \
} while(false)


void Cell:: find_intersections( const V2D &P, V2D &Q, const V2D &vmin, const V2D &vmax, const U2D &pos, Bubble *bubble )
{
    assert(bubble);
    const unit_t i  = pos.x;
    const unit_t j  = pos.y;
    const unit_t i1 = i+1;
    const unit_t j1 = j+1;
    
    int code = __OutCode(Q, vmin, vmax);
    if( INSIDE != code )
    {
        //----------------------------------------------------------------------
        // there is an intersection
        //----------------------------------------------------------------------
        if( code & TOP )
        {
            assert( Q.y > P.y );
            const Real Ix = P.x + ( vmax.y - P.y)* (Q.x - P.x) / (Q.y - P.y ) ;
            const Real Iy = vmax.y;
            __NEW_INTER(horz_seg,j1,x,X);
        }
        
        if( code & BOTTOM )
        {
            assert( Q.y < P.y );
            const Real Ix = P.x + ( vmin.y - P.y)* (Q.x - P.x) /(Q.y - P.y );
            const Real Iy = vmin.y;
            __NEW_INTER(horz_seg,j,x,X);
        }
        
        if( code & RIGHT )
        {
            assert( Q.x > P.x );
            const Real Iy = P.y + ( vmax.x - P.x )  * (Q.y - P.y ) / (Q.x - P.x );
            const Real Ix = vmax.x;
            __NEW_INTER(vert_seg,i1,y,Y);
        }
        
        if( code & LEFT )
        {
            assert( Q.x < P.x );
            const Real Iy = P.y + ( vmin.x - P.x )  * (Q.y - P.y ) / (Q.x - P.x );
            const Real Ix = vmin.x;
            __NEW_INTER(vert_seg,i,y,Y);
        }
        
        
           
        
    }
    
}


void Cell:: locate_point( Point &p )
{
    //==========================================================================
    // first pass: locate in which grid we are
    //==========================================================================
    
    //--------------------------------------------------------------------------
    // by simple bissection
    //--------------------------------------------------------------------------
    const V2D    P = p.vertex;
    assert(P.x>=0);
    assert(P.x<=Length.x);
    assert(P.y>=-Length.y/2);
    assert(P.y<=Length.y/2);
    const unit_t i = p.pos.x = __locate_point(P.x, X);
    const unit_t j = p.pos.y = __locate_point(P.y, Y);
    
    
    fprintf( stderr, " (%g,%g) <= (%g,%g) <= (%g,%g)\n", X[i],Y[j], P.x, P.y, X[i+1], Y[j+1]);
    
    //--------------------------------------------------------------------------
    // compute coefficient of bilinear interpolations
    //--------------------------------------------------------------------------
    p.w.x = (P.x - X[i])/dX[i];
    p.w.y = (P.y - Y[j])/dY[j];
    
    const V2D vmin( X[i],   Y[j] );
    const V2D vmax( X[i+1], Y[j+1]);
    
    
    assert( INSIDE == __OutCode(P,vmin,vmax) );
#if 0
    if( INSIDE != __OutCode(P, vmin, vmax) )
    {
        fprintf( stderr, "(%g,%g) is not inside (%g,%g) -> (%g,%g)\n", P.x, P.y, vmin.x, vmin.y, vmax.x, vmax.y);
        abort();
    }
#endif
    
    //--------------------------------------------------------------------------
    // find the intersections
    //--------------------------------------------------------------------------
    V2D Q    = P + p.r_next;
    find_intersections(P, Q, vmin, vmax, p.pos, p.bubble);
    
}

#include "yocto/comparator.hpp"

static inline int __compare_horz( const Segment *lhs, const Segment *rhs, void * )
{
    return __compare(lhs->inter->vertex.x,rhs->inter->vertex.x);
}

static inline int __compare_vert( const Segment *lhs, const Segment *rhs, void * )
{
    return __compare(lhs->inter->vertex.y,rhs->inter->vertex.y);
}

static inline int __compare_segment( const Segment *lhs, const Segment *rhs, void * )
{
    return lhs->inter->lo - rhs->inter->lo;
}

#include "yocto/core/merge-sort.hpp"

//! locate all points in all spots
void Cell:: locate_points( )
{
    //--------------------------------------------------------------------------
    //! bubble locator: set to zero
    //--------------------------------------------------------------------------
    B.ldz();
    Bvis.ldz();
    
    //--------------------------------------------------------------------------
    //! empty every segment
    //--------------------------------------------------------------------------
    for( size_t i=segments.size();i>0;--i) 
        segments[i].empty();
    
    //--------------------------------------------------------------------------
    //! empty intersections
    //--------------------------------------------------------------------------
    inter.empty();
    
    //--------------------------------------------------------------------------
    //! find all the intersections
    //--------------------------------------------------------------------------
    //fprintf(stderr,"....find all intersections\n");
    for( Bubble *b = bubbles.first(); b; b=b->next )
    {
        //fprintf( stderr, "#points = %lu, #spots=%lu\n", b->size, b->spots.size);
        //-- clean markers
        b->markers.empty();
        for( Spot *s = b->spots.head; s; s=s->next )
        {
            locate_point( * (s->point) );
        }
    }
    
    //--------------------------------------------------------------------------
    //! lower indices are already computed
    //--------------------------------------------------------------------------
    for( size_t i=segments.size();i>0;--i)
    {
        core::merging<Segment>::sort<core::list_of>( segments[i], __compare_segment,0);
    }
    
    //--------------------------------------------------------------------------
    //! scan horizontal segment to determine bubble segmentation
    //--------------------------------------------------------------------------
    for( unit_t j=upper.y;j>=lower.y;--j)
    {
        //----------------------------------------------------------------------
        // take the segment j
        //----------------------------------------------------------------------
        const Segment::List &Sj = horz_seg[j];
#if 0
        fprintf( stderr, "segment[%lu]@y=%8.2f: #%lu",j,Y[j],Sj.size);
        for( const Segment *seg = Sj.head; seg; seg=seg->next)
        {
            Intersection *I = seg->inter;
            assert(I);
            fprintf( stderr, " (%.3f,%.3f)/[%.3f:%3f]", I->vertex.x, I->vertex.y, X[I->lo], X[I->up]);
        }
        fprintf(stderr, "\n");
#endif
        
        //----------------------------------------------------------------------
        //  fill B[j] and Bvis[j]
        //----------------------------------------------------------------------
        ArrayInt1D &Bj  = B[j];
        Array1D    &Vj  = Bvis[j];
        unit_t   i      = Bj.lower;
        Segment *s      = Sj.head;
        bool     inside = false;
        while(s)
        {
            size_t count = 1;
            while( s->next && s->next->inter->lo == s->inter->lo )
            {
                ++count;
                s=s->next;
            }
            assert(s);
            
            //------------------------------------------------------------------
            // forward i->lo
            //------------------------------------------------------------------
            if( inside )
            {
                Bubble      *bubble    = s->inter->bubble;     assert(bubble!=NULL);
                const unit_t bubble_id = s->inter->bubble->id;
                while( i <= s->inter->lo )
                {
                    Vj[i] = Bj[i] = bubble_id;
                    const U2D marker_pos(i,j);
                    bubble->markers.append()->pos = marker_pos;
                    ++i;
                }
            }
            
            //------------------------------------------------------------------
            // update status
            //------------------------------------------------------------------
            if( (count&1) )
                inside = !inside;
            i = s->inter->up;
            s = s->next;
        }
        
    }
    // TODO: check afterwards
    
    
    
    
    
}

