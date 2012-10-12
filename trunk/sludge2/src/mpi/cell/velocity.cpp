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


void Cell:: compute_spot_velocities()
{
#if 0
    for( Bubble *bubble = bubbles.first();bubble;bubble=bubble->next)
    {
        bubble->save_vtk( vformat("b%u.vtk",bubble->id) );
        bubble->save_vtk_shell( vformat("b%u-shell.vtk",bubble->id) );
    }
#endif
   
#if 0
    for( Bubble *bubble = bubbles.first();bubble;bubble=bubble->next)
    {
        bubble->save_dat(vformat("b%u.dat",bubble->id));
    }
   
    save_B( "inside.dat" );
    {
        ios::ocstream fp("h0.dat",false);
    }
    
    {
        ios::ocstream fp("h1.dat",false);
    }
    
    
    
    ios::ocstream fp2("s.dat",false);
#endif
    
    for( Bubble *bubble = bubbles.first();bubble;bubble=bubble->next)
    {
        
        for( Spot *spot = bubble->spots.head;spot;spot=spot->next)
        {
            compute_spot_velocity(spot);
        }
        bubble->save_vtk( vformat("bubble%u.vtk", bubble->id) );
        bubble->save_vtk_n( vformat("curv%u.vtk", bubble->id) );
        bubble->save_vtk_g( vformat("bgradp%u.vtk", bubble->id) );
    }
}

// neighbors

struct neighbor
{
    Real   score;
    unit_t i;  
    unit_t j;
    friend inline bool operator<( const neighbor &lhs, const neighbor &rhs )
    {
        return rhs.score < lhs.score;
    }
};

#include "yocto/code/gsort.hpp"


void Cell:: compute_spot_velocity( Spot *spot )
{
    Tracer       *tracer = spot->handle; assert(tracer->bubble);
    const Bubble *bubble = tracer->bubble;
    
    //--------------------------------------------------------------------------
    // evaluate pressure on the tracer
    //--------------------------------------------------------------------------
    const Real    P0     = bubble->pressure -  bubble->gam * tracer->curvature;
    
    //--------------------------------------------------------------------------
    // build the probe
    //--------------------------------------------------------------------------
    const Vertex v0    = tracer->vertex;              // original point
    const Vertex vec   = -tracer->n;                  // go outside
    const Real   dx    = vec.x * delta.x;
    const Real   dy    = vec.y * delta.y;
    const Real   len   = Sqrt(dx*dx+dy*dy)*0.5;       // spacing
    const Vertex step  = len * vec;
    Vertex       probe = v0 + step;                   // starting point
    
    
    
#if 0
    {
        ios::ocstream fp("s.dat",true);
        fp("%g %g\n", probe.x, probe.y);
    }
#endif
    
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
                n.score = pg * vec; // dot product of search direction * (probe->grid)
                n.i     = I;        // where to take the pressure
                n.j     = J;        // where to take the pressure
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
    
    //--------------------------------------------------------------------------
    // Solve the pressure gradient
    //--------------------------------------------------------------------------
    const unit_t i1 = neigh[0].i;
    const unit_t j1 = neigh[0].j;
    const unit_t i2 = neigh[1].i;
    const unit_t j2 = neigh[1].j;
    //fprintf(stderr,"n0.i=%ld,n0.j=%ld,n1.i=%ld,n2.i=%ld\n",i1,j1,i2,j2);
    const Vertex v1( X[i1], Y[j1] );
    const Vertex v2( X[i2], Y[j2] );
    const Real   P1 = P[j1][i1];
    const Real   P2 = P[j2][i2];
    const Real   mA  = v1.x - v0.x;
    const Real   mB  = v1.y - v0.y;
    const Real   mC  = v2.x - v0.x;
    const Real   mD  = v2.y - v0.y;
    const Real   det = mA * mD - mB * mC;
    //fprintf(stderr,"mA=%g,mB=%g,mC=%g,mD=%g\n",mA,mB,mC,mD);
    if( Fabs(det) <= 0 )
    {
        fprintf( stderr, "invalid points for tracer@(%g,%g)\n", v0.x, v0.y);
    }
    const Real   dP1 = P1-P0;
    const Real   dP2 = P2-P0;
    const Real   gx  = ( mD * dP1 - mB * dP2)/det;
    const Real   gy  = (-mC * dP1 + mA * dP2)/det;
    spot->gradP.x = gx;
    spot->gradP.y = gy;
    spot->U = gradP_to_U( spot->gradP );
    
    
}

