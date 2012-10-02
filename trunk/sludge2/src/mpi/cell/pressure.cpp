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
        if( Pleave[j][im].x != 1)
        {
            fprintf(stderr,"Invalid Pleave[%d][%d].x=%g\n", j, im, Pleave[j][im].x );
            abort();
        }
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
        if( Penter[j][ip].x != 1)
        {
            fprintf(stderr,"Invalid Penter[%d][%d].x=%g\n", j, ip, Penter[j][ip].x );
            abort();
        }
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
        if( Pleave[jm][i].y != 1)
        {
            fprintf(stderr,"P_lower::Invalid Pleave[%d][%d].y=%g\n", jm, i, Pleave[jm][i].y );
            abort();
        }
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
        if( Penter[jp][i].y != 1)
        {
            fprintf(stderr,"P_upper::Invalid Penter[%d][%d].y=%g\n", jp, i, Penter[jp][i].y );
            abort();
        }
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
    ios::ocstream fp("py-eff.dat",false);
    
    for( unit_t i=X.lower;i<=X.upper;++i)
    {
        fp("@i=%d\n", i);
        for( unit_t j=Y.lower;j<=Y.upper;++j)
        {
            if( Penter[j][i].y > 0 )
                fp("Penter[%d][%d].y=%g\n", j, i, Penter[j][i].y);
            if( Pleave[j][i].y > 0 )
                fp("Pleave[%d][%d].y=%g\n", j, i, Pleave[j][i].y);
        }
    }
    
    
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
#if 1
                        const Real p_left  = P_left( j,i);
                        const Real p_right = P_right(j,i);
                        const Real p_lower = P_lower(j,i);
                        const Real p_upper = P_upper(j,i);
                        const Real residue =
                        inv_delsq.x * ( p_right + mid + p_left ) +
                        inv_delsq.y * ( p_upper + mid + p_lower);
#else
                        const Real residue =
                        inv_delsq.x * ( P[j][i+1] + mid + P[j][i-1] ) +
                        inv_delsq.y * ( P[j+1][i] + mid + P[j-1][i] );
#endif
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