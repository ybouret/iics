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
    Real sum_x   = 0;
    Real sum_y   = 0;
    Real sum_xx  = 0;
    Real sum_xy  = 0;
    Real sum_yy  = 0;
    Real sum_xP  = 0;
    Real sum_yP  = 0;
    J->num = 0;
    for( unit_t jj=0; jj <=1; ++jj )
    {
        const unit_t     j = klo.y + jj;
        for( unit_t ii=0; ii <=1; ++ii )
        {
            const unit_t i = klo.x + ii;
            if(B[j][i] <= 0 )
            {
                JAround &Q = J->reg[J->num++];
                const Real pressure = ( Q.P = P[j][i]);
                const Real x        = (Q.X=X[i]) - v0.x;
                const Real y        = (Q.Y=Y[j]) - v0.y;
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
    assert( J->num == Bulk[klo] );
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

#include "yocto/code/utils.hpp"

struct update_context
{
    Real sum_xdP;
    Real sum_ydP;
    Real sum_xx;
    Real sum_yy;
    Real sum_xy;
};

void update_eval( const Junction *J, const Tracer *tracer,  update_context &ctx)
{
    const Real x0 = tracer->vertex.x;
    const Real y0 = tracer->vertex.y;
    const Real P0 = tracer->pressure;
    
    for( size_t i=0; i < J->num; ++i )
    {
        const JAround & j = J->reg[i];
        const Real x  = j.X - x0;
        const Real y  = j.Y - y0;
        const Real dP = j.P - P0;
        
        ctx.sum_xdP += x * dP;
        ctx.sum_ydP += y * dP;
        ctx.sum_xx  += x*x;
        ctx.sum_xy  += x*y;
        ctx.sum_yy  += y*y;
    }
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
    const Vertex  A  = jprev->vertex;
    const Vertex  Q  = tracer->vertex;
    const Vertex  B  = jnext->vertex;
    const Vertex  AQ(A,Q);
    const Vertex  QB(Q,B);
    const Real    aq = AQ.norm();
    const Real    qb = QB.norm();
    const Real    ab = aq+qb;
    const Real    wA = clamp<Real>(0,qb/ab,1);
    const Real    wB = 1.0  - wA;
    const Real    wA2 = wA * wA;
    const Real    wB2 = wB * wB;
    
    //std::cerr << "wA=" << wA << ", wB=" << wB << std::endl;
    
    update_context ctxA,ctxB;
    memset( &ctxA, 0, sizeof(update_context));
    memset( &ctxB, 0, sizeof(update_context));
    
    update_eval(jprev, tracer, ctxA);
    update_eval(jnext, tracer, ctxB);
    
    update_context ctx;
    for( size_t i=0; i < sizeof(ctx)/sizeof(Real); ++i )
    {
        *(( (Real *)&ctx ) + i) =
        *(( (Real *)&ctxA ) + i) * wA2 +
        *(( (Real *)&ctxB ) + i) * wB2;
    }
    
    
    const Real mA = ctx.sum_xx;
    const Real mB = ctx.sum_xy;
    const Real mC = ctx.sum_xy;
    const Real mD = ctx.sum_yy;
    
    const Real vX = ctx.sum_xdP;
    const Real vY = ctx.sum_ydP;
    
    const Real det = mA * mD - mB * mC;
    
    //--------------------------------------------------------------------------
    // copy local gradient component
    //--------------------------------------------------------------------------
    const Vertex g(( mD * vX - mB * vY) / det,
                   (-mC * vX + mA * vY) / det );
    

    //--------------------------------------------------------------------------
    // keep orthonormal projection
    //--------------------------------------------------------------------------
    spot->gn = tracer->n * g;
    
    
    
    
    //--------------------------------------------------------------------------
    // gradient construction
    //--------------------------------------------------------------------------
    spot->gradP = tracer->gt * tracer->t + spot->gn * tracer->n;
    
    //--------------------------------------------------------------------------
    // velocity estimation
    //--------------------------------------------------------------------------
    spot->U     = gradP_to_U(spot->gradP);
}

