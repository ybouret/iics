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
    
    
    //segmenter.save_vtk_n( "jn.vtk");
    //segmenter.save("j.dat");
    //{ ios::ocstream fp("probe.dat",false); }
    
    for( Bubble *bubble = bubbles.first();bubble;bubble=bubble->next)
    {
        
        for( Spot *spot = bubble->spots.head;spot;spot=spot->next)
        {
            compute_spot_velocity(spot);
        }
        // bubble->save_vtk( vformat("bubble%u.vtk", bubble->id) );
        //bubble->save_vtk_gt( vformat("bgt%u.vtk", bubble->id) );
        //bubble->save_vtk_gn( vformat("bgn%u.vtk", bubble->id) );
        //bubble->save_vtk_n( vformat("curv%u.vtk", bubble->id) );
        //bubble->save_vtk_u( vformat("bu%u.vtk", bubble->id) );
        
    }
    
    //segmenter.save_vtk_gn("jgn.vtk");
    //segmenter.save_vtk_gradP( "jg.vtk" );
    
}


void Cell:: compute_junction_gn( ConstJunctionPtr J )
{
    
    //--------------------------------------------------------------------------
    // compute it only once
    //--------------------------------------------------------------------------
    if( J->visited )
        return;
    
    //--------------------------------------------------------------------------
    // locate a probe
    //--------------------------------------------------------------------------
    const Vertex v0    = J->vertex;              // original point
    const Vertex h     = -J->n;                  // go outside
    const Real   dx    = h.x * delta.x;          // ellipsoidal dx
    const Real   dy    = h.y * delta.y;          // ellipsoidal dy
    const Real   len   = Sqrt(dx*dx+dy*dy)*0.5;  // probe spacing
    const Vertex step  = len * h;                // probe step
    Vertex       probe = v0 + step;              // starting probe
    
    //{ ios::ocstream fp("probe.dat",true); fp("%g %g\n", probe.x, probe.y); }
    
    Coord klo;
    Coord khi;
    segmenter.locate_vertex(probe, klo, khi);
    if( Bulk[klo] < 2 )
    {
        fprintf( stderr, "junction @(%g,%g): bulk=%g for probe@(%g,%g)\n", v0.x, v0.y, Bulk[klo], probe.x, probe.y );
        fflush(stderr);
        abort();
    }
    
    //--------------------------------------------------------------------------
    // Least Square Fitting
    //--------------------------------------------------------------------------
    size_t count = 0;
    Real sum_x   = 0;
    Real sum_y   = 0;
    Real sum_xx  = 0;
    Real sum_xy  = 0;
    Real sum_yy  = 0;
    Real sum_xP  = 0;
    Real sum_yP  = 0;
    for( unit_t jj=0; jj <=1; ++jj )
    {
        const unit_t     j = klo.y + jj;
        for( unit_t ii=0; ii <=1; ++ii )
        {
            const unit_t i = klo.x + ii;
            if(B[j][i] <= 0 )
            {
                ++count;
                const Real pressure = P[j][i];
                const Real x        = X[i] - v0.x;
                const Real y        = Y[j] - v0.y;
                sum_x  += x;
                sum_y  += y;
                sum_xx += x*x;
                sum_xy += x*y;
                sum_yy += y*y;
                sum_xP += pressure * x;
                sum_yP += pressure * y;
            }
        }
    }
    assert( count == Bulk[klo] );
    const Real mA = sum_xx;
    const Real mB = sum_xy;
    const Real mC = sum_xy;
    const Real mD = sum_yy;
    
    const Real P0 = J->pressure;
    const Real vX = sum_xP - P0 * sum_x;
    const Real vY = sum_yP - P0 * sum_y;
    
    const Real det = mA * mD - mB * mC;
    
    //--------------------------------------------------------------------------
    // copy local gradient component
    //--------------------------------------------------------------------------
    const Vertex g(( mD * vX - mB * vY) / det,
                   (-mC * vX + mA * vY) / det );
    
    J->gn = g * J->n;
    
    
    J->visited = true;
    return;
    
    
}


void Cell:: compute_spot_velocity( Spot *spot )
{
    
    //--------------------------------------------------------------------------
    // localizing junctions
    //--------------------------------------------------------------------------
    ConstJunctionPtr jprev = 0;
    ConstJunctionPtr jnext = 0;
    segmenter.find_bracketing_junctions(jprev,jnext,spot);
    
    //--------------------------------------------------------------------------
    // sanity check
    //--------------------------------------------------------------------------
#if !defined(NDEBUG)
#if JUNCTION_TAG == 1
    if( jprev->kind == Junction::Vert )
    {
        assert(jprev->klo>Y.lower);
        assert(jprev->khi<Y.upper);
    }
#endif
#endif
    
    //--------------------------------------------------------------------------
    // compute (once) their normal gradient
    //--------------------------------------------------------------------------
    compute_junction_gn(jprev);
    compute_junction_gn(jnext);
    
    //--------------------------------------------------------------------------
    // get the tracer location
    //--------------------------------------------------------------------------
    const Tracer *tracer = spot->handle;
     
    ArcPoint A( jprev->vertex,  jprev->t,  jprev->pressure,  jprev->gt,  jprev->gn);
    ArcPoint Q( tracer->vertex, tracer->t, tracer->pressure, tracer->gt, 0 );
    ArcPoint B( jprev->vertex,  jnext->t,  jnext->pressure,  jnext->gt,  jnext->gn);

    Arc arc(A,Q,B);
    arcsol(arc);
    
    spot->gn = 0;
    
    //--------------------------------------------------------------------------
    // gradient construction
    //--------------------------------------------------------------------------
    spot->gradP = tracer->gt * tracer->t + spot->gn * tracer->n;
    
    //--------------------------------------------------------------------------
    // velocity estimation
    //--------------------------------------------------------------------------
    spot->U     = gradP_to_U(spot->gradP);
}

