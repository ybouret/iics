#include "../cell.hpp"

Real Cell:: P_left( unit_t j, unit_t i) const throw()
{
    assert(i>X.lower);
    assert(j>=Y.lower);
    assert(j<=Y.upper);
    assert(B[j][i] <= 0 ); // in bulk
    const unit_t im = i-1;
    if( B[j][im] <= 0 )
    {
        // in bulk as well
        return P[j][im];
    }
    else
    {
        // leaving a bubble along x
        return Pleave[j][im].x;
    }
}

Real Cell:: P_right( unit_t j, unit_t i) const throw()
{
    assert(i<X.upper);
    assert(j>=Y.lower);
    assert(j<=Y.upper);
    assert(B[j][i] <= 0 ); // in bulk
    const unit_t ip = i+1;
    if( B[j][ip] <= 0 )
    {
        // in bulk as well
        return P[j][ip];
    }
    else
    {
        // entering a bubble along x
        return Penter[j][ip].x;
    }
}

Real Cell:: P_lower( unit_t j, unit_t i) const throw()
{
    assert(i>=X.lower);
    assert(i<=X.upper);
    assert(j>Y.lower);
    assert(B[j][i]<=0); // in bulk
    const unit_t jm = j-1;
    if(B[jm][i] <= 0 )
    {
        // in bulk as well
        return P[jm][i];
    }
    else
    {
        // leaving a bubble along y
        return Pleave[jm][i].y;
    }
}

Real Cell:: P_upper( unit_t j, unit_t i) const throw()
{
    assert(i>=X.lower);
    assert(i<=X.upper);
    assert(j<Y.upper);
    assert(B[j][i]<=0); // in bulk
    const unit_t jp = j+1;
    if(B[jp][i] <= 0 )
    {
        // in bulk as well
        return P[jp][i];
    }
    else
    {
        // entering a bubble along y
        return Penter[jp][i].y;
    }
}


void Cell:: compute_pressure(const mpi &MPI )
{
    //==========================================================================
    //
    // Boundary conditions: initial pressure
    //
    //==========================================================================
    
        
    //--------------------------------------------------------------------------
    // boundary conditions
    //--------------------------------------------------------------------------
    for( unit_t j=lower.y;j<=upper.y;++j)
    {
        P[j][upper.x] = 0.5;
    }
    
    //----------------------------------------------------------------------
    // sync pressure
    //----------------------------------------------------------------------
    sync1(MPI,P);
    
    
    static const size_t shift[2] = { 1, 2};
    const Real ftol = 1e-5;
    for(size_t iter=1;;++iter)
    {
        
        //----------------------------------------------------------------------
        // Red/Black
        //----------------------------------------------------------------------
        int cvg = 1;
        for(size_t c=0;c<2;++c)
        {
            compute_gradP();
            //------------------------------------------------------------------
            // core
            //------------------------------------------------------------------
            for( unit_t j=lower.y;j<=upper.y;++j)
            {
                Array1D       &P_j = P[j];
                const Array1D &B_j = B[j];
                for(unit_t i=lower.x+shift[c];i<upper.x;i += 2)
                {
                    if(B_j[i]<=0)
                    {
                        Real      &P_ji    = P_j[i];
                        const Real P0      = P_ji;
                        const Real mid     = -(P0+P0);
                        const Real p_left  = P_left(j, i);
                        const Real p_right = P_right(j,i);
                        const Real p_lower = P_lower(j,i);
                        const Real p_upper = P_upper(j,i);
                        const Real residue =
                        inv_delsq.x * ( p_right + mid + p_left ) +
                        inv_delsq.y * ( p_upper + mid + p_lower);
                        const Real delta_P = -residue * rb_factor;
                        P_j[i] += delta_P;
                        if( Fabs(delta_P) > ftol * Fabs(P_ji) )
                            cvg = 0;
                    }
                }
            }
            
            //------------------------------------------------------------------
            // sides
            //------------------------------------------------------------------
            for( unit_t j=lower.y;j<=upper.y;++j)
            {
                Array1D       &P_j = P[j];
                P_j[lower.x] = (4.0*P_j[lower.x+1] - P_j[lower.x+2])/3.0;
            }
            sync1(MPI,P);
        }
        
        
        int converged = 0;
        MPI.Allreduce(&cvg, &converged, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if( MPI.CommWorldSize == converged)
        {
            MPI.Printf0(stderr, "\tcomputed pressure...\n");
            break;
        }
    }
    
}
