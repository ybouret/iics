#include "cell.hpp"

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
    return ilo;
}


static const int INSIDE = 0;  // 0000
static const int LEFT   = 1;  // 0001
static const int RIGHT  = 2;  // 0010
static const int BOTTOM = 4;  // 0100
static const int TOP    = 8;  // 1000

static inline 
int __OutCode( const V2D &q, const V2D &vmin, const V2D &vmax ) throw()
{
    int ans = INSIDE;
    if( q.x < vmin.x ) 
        ans |= LEFT;
    if( q.x > vmax.x )
        ans |= RIGHT;
    if( q.y < vmin.y )
        ans |= BOTTOM;
    if( q.y > vmax.y )
        ans |= TOP;
    return ans;
}

#define __NEW_INTER(SEGMENT,INDEX,LO,UP)  do  { \
Intersection *I = inter.append(); \
I->vertex = Q;      \
I->lo     = LO;     \
I->up     = UP;     \
I->bubble = bubble; \
SEGMENT[INDEX].append()->inter = I; \
} while(false)

//==============================================================================
// Cohen-Sutherland algorithm
//==============================================================================
void Cell:: find_intersections( const V2D &P, V2D &Q, const V2D &vmin, const V2D &vmax, const U2D &pos, Bubble *bubble )
{
    const unit_t i = pos.x;
    const unit_t j = pos.y;
    const unit_t i1 = i+1;
    const unit_t j1 = j+1;
    
    int code = __OutCode(Q, vmin, vmax);
    if( INSIDE != code )
    {
        //----------------------------------------------------------------------
        // there is an intersection
        //----------------------------------------------------------------------
    MOVE_Q:
        if( code & TOP )
        {
            assert( Q.y > P.y );
            Q.x = P.x + ( vmax.y - P.y)* (Q.x - P.x) / (Q.y - P.y ) ;
            Q.y = vmax.y;
            __NEW_INTER(horz_seg,j1,i,i1);
            goto TEST_Q;
        }
        
        if( code & BOTTOM )
        {
            assert( Q.y < P.y );
            Q.x = P.x + ( vmin.y - P.y)* (Q.x - P.x) /(Q.y - P.y );
            Q.y = vmin.y;
            __NEW_INTER(horz_seg,j,i,i1);
            goto TEST_Q;
        }
        
        if( code & RIGHT )
        {
            assert( Q.x > P.x );
            Q.y = P.y + ( vmax.x - P.x )  * (Q.y - P.y ) / (Q.x - P.x );
            Q.x = vmax.x;
            __NEW_INTER(vert_seg,i1,j,j1);
            goto TEST_Q;
        }
        
        if( code & LEFT )
        {
            assert( Q.x < P.x );
            Q.y = P.y + ( vmin.x - P.x )  * (Q.y - P.y ) / (Q.x - P.x );
            Q.x = vmin.x;
            __NEW_INTER(vert_seg,i,j,j1);
            //goto TEST_Q;
        }
        
        
    TEST_Q:
        if( INSIDE != ( code = __OutCode(Q, vmin, vmax) ) )
            goto MOVE_Q;
        
        
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
    
    
    //fprintf( stderr, " (%g,%g) <= (%g,%g) <= (%g,%g)\n", X[i],Y[j], P.x, P.y, X[i+1], Y[j+1]);
    
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
    
    
    V2D Q    = P + p.r_next;
    find_intersections(P, Q, vmin, vmax, p.pos, p.bubble);
    
}


//! locate all points in all spots
void Cell:: locate_points( )
{
    //! bubble locator: set to zero
    B.ldz();
    
    //! empty every segment
    for( size_t i=segments.size();i>0;--i) 
        segments[i].empty();
    
    //! empty intersections
    inter.empty();
    
    for( Bubble *b = bubbles.first(); b; b=b->next )
    {
        fprintf( stderr, "#points = %lu, #spots=%lu\n", b->size, b->spots.size);
        for( Spot *s = b->spots.head; s; s=s->next )
        {
            locate_point( * (s->point) );
        }
    }
    
    
    
}

