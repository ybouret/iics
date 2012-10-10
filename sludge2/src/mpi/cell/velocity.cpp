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
    {
        ios::ocstream fp("h0.dat",false);
    }
    
    {
        ios::ocstream fp("h1.dat",false);
    }

    
    
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
    inline neighbor(const Vertex org,
                    const unit_t  i,
                    const unit_t  j,
                    const Array1D &X,
                    const Array1D &Y,
                    const Vertex vec) :
    v( org ),
    k(i,j),
    m( X[i], Y[j]),
    d(v,m),
    score( vec * d )
    {
    }
    
    const Vertex v; //!< this position
    const Coord  k; //!< logical position
    const Vertex m; //!< other, one of the vertices
    const Vertex d; //!< vm, vector to get it
    const Real   score;
    
    inline ~neighbor() throw() {}
    inline neighbor( const neighbor &other ) throw() :
    v( other.v ),
    m( other.m ),
    d( other.d ),
    score( other.score )
    {
        
    }
    
    // decreasing order
    friend inline
    bool operator<( const neighbor &lhs, const neighbor &rhs )
    {
        return rhs.score < lhs.score;
    }
    
private:
    YOCTO_DISABLE_ASSIGN(neighbor);
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
    // find the position of the probe
    //--------------------------------------------------------------------------
    Coord klo;
    Coord kup;
    segmenter.locate_vertex(org,klo,kup);
    
    //--------------------------------------------------------------------------
    // find neighbors vectors and score
    //--------------------------------------------------------------------------
    neighbor nreg[] =
    {
        neighbor(org,klo.x,klo.y,X,Y,vec),
        neighbor(org,klo.x,kup.y,X,Y,vec),
        neighbor(org,kup.x,kup.y,X,Y,vec),
        neighbor(org,kup.x,klo.y,X,Y,vec)
    };
    const size_t nnum = sizeof(nreg)/sizeof(nreg[0]);
    c_sort(nreg,nnum);
    
#if 1
    for( size_t i=1; i < nnum; ++i )
    {
        assert(nreg[i-1].score>=nreg[i].score);
    }
#endif
    fprintf( stderr, "score@(%g,%g): %g %g\n", org.x, org.y, nreg[0].score,nreg[1].score);
    {
        ios::ocstream fp("h0.dat", true);
        fp("%g %g\n", nreg[0].m.x, nreg[0].m.y);
    }
    
    {
        ios::ocstream fp("h1.dat", true);
        fp("%g %g\n", nreg[1].m.x, nreg[1].m.y);
    }
}
