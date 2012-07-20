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
    
    //==========================================================================
    // second pass: use p_next to determine possible intersections
    // Cohen-Sutherland algorithm
    //==========================================================================
    assert( INSIDE == __OutCode(P,vmin,vmax) );
#if 0
    if( INSIDE != __OutCode(P, vmin, vmax) )
    {
        fprintf( stderr, "(%g,%g) is not inside (%g,%g) -> (%g,%g)\n", P.x, P.y, vmin.x, vmin.y, vmax.x, vmax.y);
        abort();
    }
#endif
    
    V2D Q    = P + p.r_next;
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
            goto TEST_Q;
        }
        
        if( code & BOTTOM )
        {
            assert( Q.y < P.y );
            Q.x = P.x + ( vmin.y - P.y)* (Q.x - P.x) /(Q.y - P.y );
            Q.y = vmin.y;
            goto TEST_Q;
        }
        
        if( code & RIGHT )
        {
            assert( Q.x > P.x );
            Q.y = P.y + ( vmax.x - P.x )  * (Q.y - P.y ) / (Q.x - P.x );
            Q.x = vmax.x;
            goto TEST_Q;
        }
        
        if( code & LEFT )
        {
            assert( Q.x < P.x );
            Q.y = P.y + ( vmin.x - P.x )  * (Q.y - P.y ) / (Q.x - P.x );
            Q.x = vmin.x;
            //goto TEST_Q;
        }
        
        
    TEST_Q:
        if( INSIDE != ( code = __OutCode(Q, vmin, vmax) ) )
            goto MOVE_Q;
        
        //----------------------------------------------------------------------
        //
        // register the intersection
        //
        //----------------------------------------------------------------------
        {
            //------------------------------------------------------------------
            // compute the intersection
            //------------------------------------------------------------------
            Point  *I = inter.create();
            I->vertex = Q;
            I->bubble = p.bubble;
            if( Q.x <= vmin.x )
            {
                // Q is on vert_seg[i]
                vert_seg[i].append(I);
            }
            
            if( Q.x >= vmax.x )
            {
                // Q is on vert_seg[i+1]
                vert_seg[i+1].append(I);
            }
            
            if( Q.y <= vmin.y )
            {
                // Q is on horz_seg[j]
                horz_seg[j].append(I);
            }
            
            if( Q.y >= vmax.y )
            {
                //Q  is on horz_seg[j+1]
                horz_seg[j+1].append(I);
            }
        }
        
    }
    
    
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
        for( Spot *s = b->spots.head; s; s=s->next )
        {
            locate_point( * (s->point) );
        }
    }
    
    
    
}

