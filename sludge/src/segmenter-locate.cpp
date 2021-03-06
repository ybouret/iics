#include "segmenter.hpp"

//==============================================================================
// bissection to locate where a point is
//==============================================================================
unit_t Segmenter::locate_point(Real           x, 
                               const Array1D &a) throw()
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


void Segmenter:: process_spot(Spot *spot)
{
    assert(segments.size()>0);
    assert(spot!=NULL);
    assert(spot->handle!=NULL);
    Tracer *p = spot->handle;
    const Vertex P = p->vertex;
    assert(P.y >= Y[Y.lower] );
    assert(P.y <  Y[Y.upper] );
    assert(P.x >= X[X.lower] );
    assert(P.x <= X[X.upper] );
    
    //--------------------------------------------------------------------------
    // find position on grid
    //--------------------------------------------------------------------------
    const unit_t i  = spot->gLower.x = locate_point( P.x, X );
    const unit_t j  = spot->gLower.y = locate_point( P.y, Y );
    const unit_t i1 = spot->gUpper.x = spot->gLower.x + 1;
    const unit_t j1 = spot->gUpper.y = spot->gLower.y + 1;
    
    //--------------------------------------------------------------------------
    // find grid boundaries
    //--------------------------------------------------------------------------
    const Vertex vmin( X[i],   Y[j]  );
    const Vertex vmax( X[i1],  Y[j1] );
    
    //--------------------------------------------------------------------------
    // compute bilinear interpolation coefficients
    //--------------------------------------------------------------------------
    spot->bary.x = (P.x-vmin.x)/dX[i];
    spot->bary.y = (P.y-vmin.y)/dY[j];
    
    //--------------------------------------------------------------------------
    // find potential intersections
    //--------------------------------------------------------------------------
    {
        const Vertex Q = P + p->edge;
        find_junctions(P, Q, vmin, vmax, spot);
    }
    
    if( !p->prev->is_spot )
    {
        const Vertex Q = P - (p->prev->edge);
        find_junctions(P, Q, vmin, vmax, spot);
    }
    
}


#define __NEW_INTER(SEGMENT,INDEX,COORD,ARR) \
do  {                                \
Junction *I = junctions.append();    \
I->vertex.x = Ix;                    \
I->vertex.y = Iy;                    \
I->bubble   = bubble;                \
SEGMENT[INDEX].append()->handle = I; \
I->up = 1 + (I->lo = locate_point( I->vertex.COORD, ARR )); \
} while(false)


void Segmenter:: find_junctions(const Vertex &P, 
                                const Vertex &Q, 
                                const Vertex &vmin, 
                                const Vertex &vmax, 
                                Spot         *spot )
{
    Tracer *p = spot->handle;
    const unit_t i  = spot->gLower.x;
    const unit_t j  = spot->gLower.y;
    const unit_t i1 = spot->gUpper.x;
    const unit_t j1 = spot->gUpper.y;
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
