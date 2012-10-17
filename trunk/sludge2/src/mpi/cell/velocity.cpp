#include "../cell.hpp"

Vertex Cell:: gradP_to_U( const Vertex &g ) const
{
    return -g;
}

void Cell:: compute_bulk_velocities()
{
    
    for(unit_t j=lower.y;j<=upper.y;++j)
    {
        VertexArray1D       &u_j = U[j];
        const VertexArray1D &g_j = gradP[j];
        for(unit_t i=lower.x;i<=upper.x;++i)
        {
            Vertex       &u = u_j[i];
            const Vertex &g = g_j[i];
            
            u = gradP_to_U(g);
        }
    }
    
}




void Cell:: compute_spots_velocity()
{
    
    
    
    
    segmenter.save_vtk_gt("gt.vtk");
    segmenter.save_vtk_n( "jn.vtk");
    
    for( Bubble *bubble = bubbles.first();bubble;bubble=bubble->next)
    {
        
        for( Spot *spot = bubble->spots.head;spot;spot=spot->next)
        {
            compute_spot_velocity(spot);
        }
        bubble->save_vtk( vformat("bubble%u.vtk", bubble->id) );
#if 0
        bubble->save_vtk_n( vformat("curv%u.vtk", bubble->id) );
        bubble->save_vtk_g( vformat("bgradp%u.vtk", bubble->id) );
#endif
    }
}

// neighbors

struct neighbor
{
    Real   score;
    unit_t i;
    unit_t j;
    Real   dx; //!< to tracer
    Real   dy; //!< to tracer
    Real   P;  //!< its pressure
    friend inline bool operator<( const neighbor &lhs, const neighbor &rhs )
    {
        return rhs.score < lhs.score;
    }
};

#include "yocto/code/gsort.hpp"

void Cell:: compute_junction_gn( Junction *J )
{
    assert(J->kind==Junction::Vert || J->kind == Junction::Horz);
}

void Cell:: compute_junctions_gn()
{
    for( unit_t j=Y.lower;j<=Y.upper;++j)
    {
        Junctions &junctions = segmenter.Horz(j);
        for( Junction *J = junctions.head;J;J=J->next)
        {
            compute_junction_gn(J);
        }
    }
    
    for(unit_t i=X.lower;i<=X.upper;++i)
    {
        Junctions &junctions = segmenter.Vert(i);
        for( Junction *J = junctions.head;J;J=J->next)
        {
            compute_junction_gn(J);
        }
    }
    
    {
        Junctions &junctions = segmenter.duplicated();
        for( Junction *J = junctions.head;J;J=J->next)
        {
            compute_junction_gn(J);
        }
    }
}


void Cell:: compute_spot_velocity( Spot *spot )
{
    
    //--------------------------------------------------------------------------
    // compute the junctions ortho gradP
    //--------------------------------------------------------------------------
    compute_junctions_gn();
    
    //--------------------------------------------------------------------------
    // localizing junctions
    //--------------------------------------------------------------------------
    ConstJunctionPtr jprev = 0;
    ConstJunctionPtr jnext = 0;
    segmenter.find_bracketing_junctions(jprev,jnext,spot);
    
    
#if 0
    Tracer       *tracer = spot->handle; assert(tracer->bubble);
    const Bubble *bubble = tracer->bubble;
    
    //--------------------------------------------------------------------------
    // evaluate pressure on the tracer
    //--------------------------------------------------------------------------
    const Real    P0     = bubble->pressure -  bubble->gam * tracer->curvature;
    
    //--------------------------------------------------------------------------
    // build the probe : ellipsoidal direction
    //--------------------------------------------------------------------------
    const Vertex v0    = tracer->vertex;              // original point
    const Vertex h     = -tracer->n;                  // go outside
    const Real   dx    = h.x * delta.x;
    const Real   dy    = h.y * delta.y;
    const Real   len   = Sqrt(dx*dx+dy*dy)*0.5;       // spacing
    const Vertex step  = len * h;
    Vertex       probe = v0 + step;                   // starting point
    
    
    
    //--------------------------------------------------------------------------
    // find the right position
    //--------------------------------------------------------------------------
    Coord klo;
    Coord kup;
    for(;;)
    {
        //----------------------------------------------------------------------
        // check the probe position
        //----------------------------------------------------------------------
        if(probe.x<=0) probe.x = 0;
        
        segmenter.locate_vertex(probe,klo,kup);
        if( Bulk[klo] < 2 )
        {
            //fprintf( stderr, "tracer @(%g,%g): invalid probe (bulk lo@(%g,%g)=%g)\n", v0.x, v0.y, X[klo.x], Y[klo.y], Bulk[klo.y][klo.x]);
            probe += step;
            continue;
        }
        assert(Bulk[klo]>=2);
        break;
    }
    
    if(false)
    {
        ios::ocstream fp("s.dat",true);
        fp("%g %g\n", probe.x, probe.y);
    }
    
    //--------------------------------------------------------------------------
    // find the right neigbhors
    //--------------------------------------------------------------------------
    neighbor  neigh[4];
    size_t    count = 0;
    for( unit_t j=0;j<2;++j)
    {
        const unit_t J = klo.y+j;
        const Real   y = Y[J];
        for(unit_t i=0;i<2;++i)
        {
            const unit_t I=klo.x+i;
            if( B[J][I] <= 0)
            {
                const Real   x=X[I];
                const Vertex g(x,y);           // grid vertex
                const Vertex pg(probe,g);      // probe -> grid
                assert(count<4);               // sanity check
                neighbor &n = neigh[count];
                n.score     = pg * h;   // dot product of search direction * (probe->grid)
                n.i         = I;        // where to take the pressure
                n.j         = J;        // where to take the pressure
                n.P         = P[J][I];
                n.dx        = X[I] - v0.x;
                n.dy        = Y[J] - v0.y;
                ++count;
            }
        }
    }
    assert( Bulk[klo] == count );
    assert(count>=2);
    c_sort(neigh, count);
#if !defined(NDEBUG)
    for( size_t i=1; i < count;++i) {  assert(neigh[i-1].score >= neigh[i].score); }
#endif
    
    Real sum_xP = 0;
    Real sum_yP = 0;
    Real sum_x  = 0;
    Real sum_y  = 0;
    Real sum_xx = 0;
    Real sum_yy = 0;
    Real sum_xy = 0;
    for( size_t i=0; i < count; ++i )
    {
        const neighbor &n = neigh[i];
        sum_xP += n.dx * n.P;
        sum_yP += n.dy * n.P;
        sum_x  += n.dx;
        sum_y  += n.dy;
        sum_xx += n.dx * n.dx;
        sum_yy += n.dy * n.dy;
        sum_xy += n.dx * n.dy;
    }
    
    const Real det = sum_xx * sum_yy - sum_xy * sum_xy;
    if( Fabs(det) <= 0 )
    {
        fprintf( stderr, "invalid points for tracer@(%g,%g)\n", v0.x, v0.y);
    }
    
    const Real Bx = sum_xP - sum_x * P0;
    const Real By = sum_yP - sum_y * P0;
    
    const Real gx = (  sum_yy * Bx - sum_xy * By)/det;
    const Real gy = ( -sum_xy * Bx + sum_xx * By)/det;
    const Real   dPdh = gx*h.x + gy*h.y;
    
    spot->gradP = dPdh * h;
    spot->U.ldz();
#endif
    
#if 0
    //--------------------------------------------------------------------------
    // now we have at least two neighbors, ranked by their dot products
    //--------------------------------------------------------------------------
    const unit_t i1 = neigh[0].i;
    const unit_t j1 = neigh[0].j;
    const unit_t i2 = neigh[1].i;
    const unit_t j2 = neigh[1].j;
    const Real   x1 = X[i1];
    const Real   x2 = X[i2];
    const Real   y1 = Y[j1];
    const Real   y2 = Y[j2];
    
    //--------------------------------------------------------------------------
    // estimate the normal gradient pressure
    //--------------------------------------------------------------------------
    const Real   dP1 = P[j1][i1] - P0;
    const Real   dP2 = P[j2][i2] - P0;
    
    const Real   mA  = x1 - v0.x;
    const Real   mB  = y1 - v0.y;
    const Real   mC  = x2 - v0.x;
    const Real   mD  = y2 - v0.y;
    const Real   det = mA * mD - mB * mC;
    //fprintf(stderr,"mA=%g,mB=%g,mC=%g,mD=%g\n",mA,mB,mC,mD);
    if( Fabs(det) <= 0 )
    {
        fprintf( stderr, "invalid points for tracer@(%g,%g)\n", v0.x, v0.y);
    }
    
    const Real   gx  = ( mD * dP1 - mB * dP2)/det;
    const Real   gy  = (-mC * dP1 + mA * dP2)/det;
    const Real   dPdh = gx*h.x + gy*h.y;
    
    spot->gradP = dPdh * h;
    spot->U.ldz();
#endif
    
#if 0
    //--------------------------------------------------------------------------
    // find the intersection with the edge
    //--------------------------------------------------------------------------
    const Real dx12 = x2-x1;
    const Real dy12 = y2-y1;
    
    const Real cross_h = dx12 * h.y - dy12 * h.x;
    if( Fabs(cross_h) <= 0 )
    {
        fprintf( stderr, "invalid points for tracer@(%g,%g)\n", v0.x, v0.y);
    }
    
    const Real cross_v = dx12 * (v0.y-y1) - dy12 * (v0.x-x1);
    const Real mu      = -cross_v/cross_h;
    if(mu<=0)
    {
        fprintf( stderr, "invalid mu=%g for tracer@(%g,%g)\n", mu, v0.x, v0.y);
    }
    assert(mu>0);
    
    const Vertex shell = v0 + mu * h;
    
    {
        ios::ocstream fp("s.dat",true);
        fp("%g %g\n", shell.x, shell.y);
    }
#endif
    
#if 0
    {
        ios::ocstream fp("h0.dat", true);
        fp("%g %g\n", X[ neigh[0].i ], Y[ neigh[0].j ] );
    }
    
    {
        ios::ocstream fp("h1.dat", true);
        fp("%g %g\n", X[ neigh[1].i ], Y[ neigh[1].j ] );
    }
#endif
    
}

