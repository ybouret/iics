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
    
    for( Bubble *bubble = bubbles.first();bubble;bubble=bubble->next)
    {
        
        for( Spot *spot = bubble->spots.head;spot;spot=spot->next)
        {
            compute_spot_velocity(spot);
        }
        
    }
}

// neighbors
#define TOP_LEFT  0x01
#define TOP_RIGHT 0x02
#define BOT_LEFT  0x04
#define BOT_RIGHT 0x08

void Cell:: compute_spot_velocity( Spot *spot )
{
    Tracer       *tracer = spot->handle; assert(tracer->bubble);
    const Bubble *bubble = tracer->bubble;
    
    //--------------------------------------------------------------------------
    // evaluate pressure on the tracer
    //--------------------------------------------------------------------------
    const Real    P0     = bubble->pressure -  bubble->gam * tracer->curvature;
    (void)P0;
    const Vertex v     = tracer->vertex;
    Vertex       probe = v - bubble->lam * tracer->n;
    pbc(probe);
    
    int bulk = 0;
    segmenter.locate_vertex(probe, spot->klo, spot->kup);
    const unit_t ilo = spot->klo.x;
    const unit_t iup = spot->kup.x;
    const unit_t jlo = spot->klo.y;
    const unit_t jup = spot->kup.y;
    
    if( B[jlo][ilo] <= 0)
        bulk |= BOT_LEFT;
    
    if( B[jlo][iup] <= 0)
        bulk |= BOT_RIGHT;
    
    if( B[jup][ilo] <= 0)
        bulk |= TOP_LEFT;
    
    if( B[jup][iup] <= 0)
        bulk |= TOP_RIGHT;
    
    if( !bulk )
    {
        throw exception("Can find bulk for tracer @(%g,%g)\n", v.x,v.y);
    }
    
}
