#include "workspace.hpp"
#include "yocto/code/unique.hpp"

Vertex Workspace:: gradP_to_V( const Vertex &g ) const
{
    return -g;
}


#define MAX_LOCAL_PRESSURES 12

static inline
bool __collectPressure(Workspace::LocalPressure &lp,
                       const unit_t   i,
                       const unit_t   j,
                       const Array   &P,
                       const Array   &B,
                       const Array1D &X,
                       const Array1D &Y)
{
    if((i>=X.lower) &&
       (i<=X.upper) &&
       (j>=Y.lower) &&
       (j<=Y.upper) &&
       (B[j][i]<0) )
    {
        lp.r.x = X[i];
        lp.r.y = Y[j];
        lp.P   = P[j][i];
        return true;
    }
    return false;
}

#define COLLECT_P(u,v) do { assert(n<MAX_LOCAL_PRESSURES); if( __collectPressure(lp[n], u, v, P, B, X, Y) ) ++n; } while(false)

void Workspace:: collect_pressure( const Junction *J, LocalPressure lp[], size_t &n) const
{
    assert(J);
    assert(lp);
    switch(J->root.type)
    {
        case Junction::Vert:
        {
            const unit_t i   = J->root.indx;
            const unit_t im  = i-1;
            const unit_t ip  = i+1;
            const unit_t jlo = J->lower;
            const unit_t jup = J->upper;
            
            COLLECT_P(i, jlo);
            COLLECT_P(i, jup);
            COLLECT_P(im,jlo);
            COLLECT_P(im,jup);
            COLLECT_P(ip,jlo);
            COLLECT_P(ip,jup);
        }
            break;
            
        case Junction::Horz:
        {
            const unit_t j   = J->root.indx;
            const unit_t jm  = j-1;
            const unit_t jp  = j+1;
            const unit_t ilo = J->lower;
            const unit_t iup = J->upper;
            
            COLLECT_P(ilo,j);
            COLLECT_P(ilo,jm);
            COLLECT_P(ilo,jp);
            COLLECT_P(iup,j);
            COLLECT_P(iup,jm);
            COLLECT_P(iup,jp);
        }
            break;
    }
}


void Workspace:: compute_velocities()
{
    
    //==========================================================================
    //
    // compute velocities in the bulk
    //
    //==========================================================================
    for(unit_t j=outline.lower.y;j<=outline.upper.y;++j)
    {
        for(unit_t i=outline.lower.x;i<=outline.upper.x;++i)
        {
            if(B[j][i] < 0 )
            {
                V[j][i] = gradP_to_V(gradP[j][i]);
            }
            else
                V[j][i].ldz();
        }
        
    }
    
    //==========================================================================
    //
    // compute for all the markers of all the bubbles
    //
    //==========================================================================
    LocalPressure lp[MAX_LOCAL_PRESSURES];
    junctions.save_dat( "j.dat" );

    for( Bubble *b = bubbles.head;b;b=b->next)
    {
        const Real P_in  = b->pressure;
        const Real gamma = b->gamma;
        const Real mu    = b->lambda / 2;
        b->save_dat( vformat("b%u.dat", unsigned(b->UID) ) );

        for( Marker *m = b->markers.head;m;m=m->next)
        {
            //------------------------------------------------------------------
            // Tangential, easy
            //------------------------------------------------------------------

            const Tracer *tr    = m->tracer;
            const Real    Pcurr = P_in - gamma * tr->C;
            const Real    Pnext = P_in - gamma * tr->next->C;
            const Real    Pprev = P_in - gamma * tr->prev->C;
            const Real    Pp    = Pnext-Pcurr;
            const Real    Pm    = Pprev-Pcurr;
            const Real    tp    = tr->dist;
            const Real    tm    = tr->prev->dist;
            const Real    Vp    = Pp/tp;
            const Real    Vm    = Pm/tm;
            m->gt = (tm*Vp - tp*Vm)/(tm+tp);
            
            //------------------------------------------------------------------
            // normal pressure
            //------------------------------------------------------------------
            m->gn = 0;
            junctions.bracket(*b, m);
            assert(m->jprev);
            assert(m->jnext);
            size_t np = 0;
            collect_pressure(m->jprev, lp, np);
            collect_pressure(m->jnext, lp, np);
            np = unique(lp, np, LocalPressure::CompareByVertex);
            std::cerr << "np=" << np << std::endl;
            const Vertex pos = tr->pos;
            if(np<=1)
                throw exception("Not enough neighbors for tracer @[%g %g]", pos.x, pos.y);

            // get the delta
            for(size_t i=0;i<np;++i)
            {
                lp[i].r -= pos;
                lp[i].d  = lp[i].r.norm();
            }
            
            // order by decreasing distance
            hsort(lp, np, LocalPressure::CompareByDecreasingDistance);
            std::cerr << "d="; for(size_t i=0;i<np;++i) std::cerr << " " << lp[i].d; std::cerr << " / lambda=" << b->lambda << std::endl;
            
            while(np>0 && lp[np-1].d<mu) --np;
            if(np<=1)
                throw exception("Neighbors are too close for tracer @[%g %g]", pos.x, pos.y);
            
            // take the least-square approximation
            const Real Pin    = b->pressure;
            const Real P0     = Pin - b->gamma * tr->C;
            Real       sum_xp = 0;
            Real       sum_yp = 0;
            Real       sum_x2 = 0;
            Real       sum_y2 = 0;
            Real       sum_xy = 0;
            for(size_t i=0;i<np;++i)
            {
                const LocalPressure &l = lp[i];
                const Real   x = l.r.x;
                const Real   y = l.r.y;
                const Real   p = l.P - P0;
                sum_xp += x*p;
                sum_yp += y*p;
                sum_x2 += x*x;
                sum_y2 += y*y;
                sum_xy += x*y;
            }
            const Real D = sum_x2 * sum_y2 - sum_xy * sum_xy;
            if( Fabs(D) <= numeric<Real>::minimum )
                throw exception("Singular Pressure Expression @[%g %g]\n", pos.x, pos.y);
            
            const Real alpha = ( sum_y2*sum_xp - sum_xy * sum_yp)/D;
            const Real beta  = (-sum_xy*sum_xp + sum_x2 * sum_yp)/D;
            
            
            // project onto normal
            m->gn = alpha * tr->n.x + beta * tr->n.y;
            std::cerr << "gt=" << m->gt << ", gn=" << m->gn << std::endl;
        }
        
      
    }
    
       
}

void Workspace:: save_markers( const mpi &MPI ) const
{
    const string &mpi_id = MPI.CommWorldID;
    for( const Bubble *b = bubbles.head;b;b=b->next)
    {
        const string bb_id = vformat("%u", unsigned(b->UID) ) + "-" + mpi_id + ".dat";
        b->save_dat( "b" + bb_id);
        
        {
            ios::ocstream fp( "jp" + bb_id,false);
            for( Marker *m = b->markers.head;m;m=m->next)
            {
                const Vertex org = m->tracer->pos;
                fp("%g %g\n", org.x, org.y);
                const Vertex j  = m->jprev->get();
                fp("%g %g\n\n", j.x, j.y);
            }
        }
        
        {
            ios::ocstream fp( "jn" + bb_id,false);
            for( Marker *m = b->markers.head;m;m=m->next)
            {
                const Vertex org = m->tracer->pos;
                fp("%g %g\n", org.x, org.y);
                const Vertex j  = m->jnext->get();
                fp("%g %g\n\n", j.x, j.y);
            }
        }

        
    }
}


