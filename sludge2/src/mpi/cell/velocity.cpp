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
    
    bubbles.first()->save_dat("b.dat");
    ios::ocstream fp("h.dat",false);
    ios::ocstream fp2("s.dat",false);
    
    for( Bubble *bubble = bubbles.first();bubble;bubble=bubble->next)
    {
        
        for( Spot *spot = bubble->spots.head;spot;spot=spot->next)
        {
            compute_spot_velocity(spot);
        }
        
    }
}

// neighbors
#include "yocto/code/utils.hpp"

class neighbor
{
public:
    inline neighbor( const Vertex org, const Real x, const Real y) :
    v( org ),
    m( x, y)
    {
    }
    
    const Vertex v; //!< this
    const Vertex m; //!< other
    
    inline ~neighbor() throw() {}
    
private:

};

void Cell:: compute_spot_velocity( Spot *spot )
{
    Tracer       *tracer = spot->handle; assert(tracer->bubble);
    const Bubble *bubble = tracer->bubble;
    
    //--------------------------------------------------------------------------
    // evaluate pressure on the tracer
    //--------------------------------------------------------------------------
    const Real    P0     = bubble->pressure -  bubble->gam * tracer->curvature;
    (void)P0;
    
    //--------------------------------------------------------------------------
    // build the probe
    //--------------------------------------------------------------------------
    const Vertex v0    = tracer->vertex;              // original point
    const Vertex vec   = -tracer->n;                  // go outside
    const Real   dx    = vec.x * delta.x;
    const Real   dy    = vec.y * delta.y;
    const Real   len   = Sqrt(dx*dx+dy*dy)*0.5;        // spacing
    Vertex       org   = v0 + len * vec;              // starting point
    //pbc(org);
    if(org.x<=0) org.x = 0;
    ios::ocstream fp2("s.dat", true );
    fp2("%g %g\n", org.x, org.y);
    
    //--------------------------------------------------------------------------
    // find the position or the probe
    //--------------------------------------------------------------------------
    Coord klo;
    Coord kup;
    segmenter.locate_vertex(org,klo,kup);
    
      
    
}
