#include "workspace.hpp"
#include "yocto/code/unique.hpp"

Vertex Workspace:: gradP_to_V( const Vertex &g ) const
{
    return -g;
}


#define MAX_LOCAL_PRESSURES 12

static inline
bool __collectPressure(LocalPressure      &lp,
                       const unit_t         i,
                       const unit_t         j,
                       const Array         &P,
                       const Array         &B,
                       const Array1D       &X,
                       const Array1D       &Y,
                       const VertexArray   &G)
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
        lp.g   = G[j][i];
        return true;
    }
    return false;
}

#define COLLECT_P(u,v) do { assert(n<MAX_LOCAL_PRESSURES); if( __collectPressure(lp[n], u, v, P, B, X, Y, gradP) ) ++n; } while(false)

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

#include "yocto/code/utils.hpp"

#if 0
static inline
Real __eval_pressure_at( const Vertex &probe, LocalPressure lp[], const size_t np)
{
    assert(lp!=0);
    assert(np>0);
    //--------------------------------------------------------------------------
    // first pass: initialize points
    //--------------------------------------------------------------------------
    
    Real D = 0;
    for(size_t i=0; i < np; ++i )
    {
        lp[i].dr = probe - lp[i].r;        //!< displacement to probe
        lp[i].d  = lp[i].dr.norm();
        D += lp[i].d; //!< cumulative length
    }
    
    
    Real p = 0;
    switch(np)
    {
        case 0:
            throw exception("Unexpected np=0");
            
        case 1:
            p = lp[0].P + lp[0].dr * lp[0].g;
            break;
            
        default:
            //-- second pass: add weighted pressure
            for(size_t i=0;i<np;++i)
            {
                const Real w = D - lp[i].d;
                p += w*(lp[i].P + lp[i].dr * lp[i].g);
            }
            p /= D*(np-1);
    }
    return p;
}

static inline
Real __eval_gradient_along( const Vertex &u, LocalPressure lp[], const size_t np, const Vertex &pos)
{
    assert(lp!=0);
    assert(np>0);

    //--------------------------------------------------------------------------
    // first pass: initialize points
    //--------------------------------------------------------------------------
    Real D = 0;
    for(size_t i=0; i < np; ++i )
    {
        lp[i].dr = pos - lp[i].r;  //!< displacement to tracer
        lp[i].d  = lp[i].dr.norm();
        D += lp[i].d; //!< cumulative length
    }

    //--------------------------------------------------------------------------
    //-- second pass: add weighted gradients
    //--------------------------------------------------------------------------
    Real g=0;
    
    switch(np)
    {
        case 0:
            throw exception("Unexpected np=0");
            
        case 1:
            g += lp[0].g * u;
            break;
            
        default:
            //-- second pass: add weighted pressure
            for(size_t i=0;i<np;++i)
            {
                const Real w = D - lp[i].d;
                g += w*(lp[i].g*u);
            }
            g /= D*(np-1);
    }

    return g;
    
}
#endif

namespace
{
    struct BadVertex
    {
        Vertex n_out;
        bool operator()( const LocalPressure &lp )
        {
            return (lp.dr * n_out) <= 0;
        }
    };
}

#include "yocto/code/remove-if.hpp"

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
    
#define SAVE_INFO 0
    
    const Vertex ex(1,0);
    for( Bubble *b = bubbles.head;b;b=b->next)
    {
        const Real P_in  = b->pressure;
        const Real gamma = b->gamma;
        
#if defined(SAVE_INFO) && SAVE_INFO == 1
        b->save_dat( vformat("b%u.dat", unsigned(b->UID) ) );
        
        
        ios::ocstream fp ( vformat("lp%u.dat", unsigned(b->UID)), false);
        ios::ocstream fp2( vformat("gn%u.dat", unsigned(b->UID)), false);
        //ios::ocstream fp3( vformat("pr%u.dat", unsigned(b->UID)), false);
#endif
        
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
            // remove invalid vertices
            //------------------------------------------------------------------
            Real D = 0;
            for(size_t i=0; i < np; ++i )
            {
                lp[i].dr = lp[i].r - pos;
                lp[i].d  = lp[i].dr.norm();
                D += lp[i].d;
            }
            
            BadVertex is_bad;
            is_bad.n_out = - tr->n;
            np = remove_if(lp, np, is_bad);
            if(np<=0)
                throw exception("Not enough VALID neighbors for tracer @[%g %g]", pos.x, pos.y);
            
            
            const Vertex t = tr->t;
            const Vertex n = tr->n;
            
            Real num = 0;
            Real den = 0;
            const Real alpha = m->gt;
            for(size_t i=0; i < np; ++i )
            {
                const Vertex &AQ     = lp[i].dr;
                const Real    dP     = lp[i].P - (Pcurr+alpha*(t*AQ));
                const Real    weight = np > 1 ? D - lp[i].d : 1;
                const Real    fac    =  (n*AQ);
                const Real    wfac   = weight * fac;
                num += wfac * dP;
                den += wfac * fac;
            }
            m->gn = num/den;

            
#if 0
            hsort(lp, np, LocalPressure::CompareByIncreasingDistance);
            //np = min_of<size_t>(np,2);
            
            Real x2 = 0;
            Real xy = 0;
            Real y2 = 0;
            Real xP = 0;
            Real yP = 0;
            
            for(size_t i=0; i < np; ++i )
            {
                const Vertex &dr = lp[i].dr;
                x2 += dr.x * dr.x;
                xy += dr.x * dr.y;
                y2 += dr.y * dr.y;
                xP += dr.x * ( lp[i].P - Pcurr);
                yP += dr.y * ( lp[i].P - Pcurr);
            }
            const Vertex t = tr->t;
            MM[1][1] = x2;  MM[1][2] = xy;  MM[1][3] = -t.x;
            MM[2][1] = xy;  MM[2][2] = y2;  MM[2][3] = -t.y;
            MM[3][1] = t.x; MM[3][2] = t.y; MM[3][3] = 0;
            
            if( !LU.build(MM) )
                throw exception("Invalid Neigborhood for tracer @[%g %g]", pos.x, pos.y);
            
            MR[1] = xP;
            MR[2] = yP;
            MR[3] = m->gt;
            
            
            LU.solve(MM, MR);
            
            const Real alpha = MR[1];
            const Real beta  = MR[2];
            m->gn = alpha * tr->n.x + beta * tr->n.y;
#endif
            
#if defined(SAVE_INFO) && SAVE_INFO == 1
            
            for( size_t i=0; i < np; ++i )
            {
                fp("%g %g\n",   pos.x, pos.y);
                fp("%g %g\n\n", lp[i].r.x, lp[i].r.y );
                
                //fp3("%g %g\n", pos.x, pos.y);
                //fp3("%g %g\n\n", probe.x, probe.y);
            }
            
            fp2("%g %g\n", pos.x+m->gn*tr->n.x, pos.y+m->gn*tr->n.y);
#endif
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


