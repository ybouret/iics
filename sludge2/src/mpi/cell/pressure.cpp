#include "../cell.hpp"

void Cell:: compute_pressure(const mpi &MPI )
{
    //==========================================================================
    //
    // Boundary conditions: initial pressure
    //
    //==========================================================================
    
    //--------------------------------------------------------------------------
    // pressure from bubbles
    //--------------------------------------------------------------------------
    segmenter.pressurize(P);
    
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
                        const Real residue =
                        inv_delsq.x * ( P_j[i+1]  + mid + P_j[i-1]) +
                        inv_delsq.y * ( P[j+1][i] + mid + P[j-1][i]);
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
