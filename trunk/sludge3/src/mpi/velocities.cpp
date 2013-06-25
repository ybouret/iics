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
            unit_t ja=0,jb=0;
            
            switch( J->b_pos )
            {
                case Bubble::IsAfter:
                    ja = J->lower;
                    jb = ja-1;
                    break;
                    
                case Bubble::IsBefore:
                    ja = J->upper;
                    jb = ja+1;
                    break;
                    
            default:
                throw exception("Invalid Vertical Bubble Position");
            }
            
            COLLECT_P(im, ja);
            COLLECT_P(i,  ja);
            COLLECT_P(ip, ja);
            
            COLLECT_P(im, jb);
            COLLECT_P(i,  jb);
            COLLECT_P(ip, jb);
            
        }
            break;
            
        case Junction::Horz:
        {
            const unit_t j   = J->root.indx;
            const unit_t jm  = j-1;
            const unit_t jp  = j+1;
            unit_t ia=0,ib=0;
            switch( J->b_pos )
            {
                    
                case Bubble::IsAfter:
                    ia = J->lower;
                    ib = ia-1;
                    break;
                    
                case Bubble::IsBefore:
                    ia = J->upper;
                    ib = ia+1;
                    break;
                    
                default:
                    throw exception("Invalid Horizontal Bubble Position");
            }
            
            COLLECT_P(ia,jm);
            COLLECT_P(ia,j);
            COLLECT_P(ia,jp);
            
            COLLECT_P(ib,jm);
            COLLECT_P(ib,j);
            COLLECT_P(ib,jp);
            
            
        }
            break;
    }
}

#include "yocto/code/remove-if.hpp"

namespace
{
    struct InvalidVertex
    {
        Vertex out_n;
        Real   d_min;
        
        bool operator()( const Workspace::LocalPressure &lp )
        {
            return (lp.r*out_n <= 0) || lp.d < d_min;
        }
    };
    
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
    //junctions.save_dat( "j.dat" );
    
    for( Bubble *b = bubbles.head;b;b=b->next)
    {
        const Real P_in  = b->pressure;
        const Real gamma = b->gamma;
        const Real mu    = b->lambda / 2;
        b->save_dat( vformat("b%u.dat", unsigned(b->UID) ) );
        
        
        ios::ocstream fp( vformat("lp%u.dat", unsigned(b->UID)) , false);
        ios::ocstream fp2( vformat("gn%u.dat", unsigned(b->UID)) , false);

        for( Marker *m = b->markers.head;m;m=m->next)
        {
            //------------------------------------------------------------------
            //
            // Tangential, easy
            //
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
            //
            // normal pressure
            //
            //------------------------------------------------------------------
            m->gn = 0;
            junctions.bracket(*b, m);
            assert(m->jprev);
            assert(m->jnext);
            
            //------------------------------------------------------------------
            // collect all the possible local pressures
            //------------------------------------------------------------------
            size_t np = 0;
            collect_pressure(m->jprev, lp, np);
            collect_pressure(m->jnext, lp, np);
            
            //------------------------------------------------------------------
            // remove doublons
            //------------------------------------------------------------------
            np = unique(lp, np, LocalPressure::CompareByVertex);
            
            const Vertex pos = tr->pos;
            if(np<=0)
                throw exception("Not enough neighbors for tracer @[%g %g]", pos.x, pos.y);
            
            //------------------------------------------------------------------
            // get the delta, compute the distance to tracer
            //------------------------------------------------------------------
            for(size_t i=0;i<np;++i)
            {
                lp[i].r -= pos;
                lp[i].d  = lp[i].r.norm();
            }
            
            InvalidVertex chk;
            chk.out_n = -tr->n;
            chk.d_min = mu;
            
            np = remove_if(lp,np,chk);
            
            if(np<=0)
                throw exception("Not enough valid neighbors for tracer @[%g %g]", pos.x, pos.y);
            
            
            
            const Real alpha = m->gt;
            
            Real         weight  = 0;
            Real         residue = 0;
            const Vertex t = tr->t;
            const Vertex n = tr->n;
            const Real   P0 = b->pressure - tr->C * b->gamma;
            for( size_t i=0; i < np; ++i )
            {
                const LocalPressure &l = lp[i];
                const Vertex        &r = l.r;
                const Real           coef = r*n;
                weight  += coef*coef;
                residue += coef*(l.P - (P0+alpha*(r*t)));
                fp("%g %g\n", pos.x, pos.y);
                fp("%g %g\n\n", pos.x + r.x, pos.y + r.y );
            }
            m->gn = residue / weight;
            //std::cerr << "np=" << np << ", gt=" << m->gt << ", gn=" << m->gn << std::endl;
            fp2("%g %g\n", pos.x, pos.y);
            fp2("%g %g\n\n", pos.x+m->gn*tr->n.x, pos.y+m->gn*tr->n.y);
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
            ios::ocstream fp("m" + bb_id, false);
            for( Marker *m = b->markers.head;m;m=m->next)
            {
                const Vertex org = m->tracer->pos;
                fp("%g %g\n", org.x, org.y);
            }
        }
        
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


