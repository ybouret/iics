#include "../cell.hpp"

static inline
const Junction *find_first_at( const unit_t klo, const Junctions &junctions)
{
    const Junction *J = junctions.head;
    
    //-- find first corresponding junction
    while( J )
    {
        if(J->klo==klo)
            break;
        J=J->next;
    }
    if( !J )
        throw exception("can't find junction for k=%d !", int(klo));
    
    //-- done
    return J;
}

static inline
const Junction *find_last_at( const unit_t klo, const Junctions &junctions)
{
    const Junction *J = find_first_at(klo,junctions);
    
    //-- forward
    while(J->next&&J->next->klo==klo)
        J=J->next;
    
    //-- done
    return J;
}


static inline
Real derivatives( const Real z[], const Real f[] )
{
    assert(z);
    assert(f);
    
    return 0;
}


void Cell:: compute_gradP()
{
    
    //--------------------------------------------------------------------------
    // initialize: ok for walls and inside bubbles
    //--------------------------------------------------------------------------
    gradP.ldz();
    
    
    //--------------------------------------------------------------------------
    // vertical gradient
    //--------------------------------------------------------------------------
    Real z[3];
    Real f[3];
    
    for( unit_t i=lower.x+1;i<upper.x;++i)
    {
        const Junctions &junctions = segmenter.Vert(i);
        for( unit_t j=lower.y;j<=upper.y;++j)
        {
            if(B[j][i]>0)
                continue; // in a bubble
            
            assert(B[j][i]<=0);
            //------------------------------------------------------------------
            // record the central point
            //------------------------------------------------------------------
            z[1] = Y[j];
            f[1] = P[j][i];
            
            //------------------------------------------------------------------
            // record the previous point
            //------------------------------------------------------------------
            const unit_t jm = j-1;
            if( B[jm][i] <= 0 )
            {
                z[0] = Y[jm];
                f[0] = P[jm][i];
            }
            else
            {
                const Junction *J = find_last_at(jm,junctions);
                z[0] = J->vertex.y;
                f[0] = J->bubble->pressure; // +gamma*J->curvature
            }
            //------------------------------------------------------------------
            // record the next point
            //------------------------------------------------------------------
            const unit_t jp = j+1;
            if(B[jp][i]>=0)
            {
                z[2] = Y[jp];
                f[2] = P[jp][i];
            }
            else
            {
                const Junction *J = find_first_at(j,junctions);
                z[2] = J->vertex.y;
                f[2] = J->bubble->pressure; // +gamma*J->curvature
            }
            
            gradP[j][i].y = derivatives(z,f);
            
        }
    }
    
    
}
