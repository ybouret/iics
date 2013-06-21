#include "workspace.hpp"
#include "yocto/mpi/ops.hpp"

void Workspace:: compute_laplacian( )
{
    DeltaP.ldz();
    for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
    {
        const unit_t jm = j-1;
        const unit_t jp = j+1;
        for( unit_t i=bulk_imin; i <= bulk_imax; ++i )
        {
            if(B[j][i]<0)
            {
                const Real P_center  = P[j][i];
                const Real P_left    = L2[j][i-1].x;
                const Real P_right   = E2[j][i+1].x;
                const Real P_bottom  = L2[jm][i].y;
                const Real P_top     = E2[jp][i].y;
                const Real mid       = -(P_center+P_center);
                const Real Laplacian = (P_left+mid+P_right) * order2fac.x + (P_bottom+mid+P_top) * order2fac.y;
                DeltaP[j][i] = Laplacian;
            }
            
        }
        
    }
    
}

void Workspace:: compute_pressure( const mpi &MPI )
{
    const int target = MPI.CommWorldSize;
    for(;;)
    {
        const int cvg = update_pressure(MPI, Red) & update_pressure(MPI, Black);
        if( target == mpi_ops::sum_of(MPI, cvg, MPI_COMM_WORLD) )
            break;
    }
    
    MPI.Printf0(stderr,"\t\tpressurize bubbles...\n");
    pressurize_bubbles();
    
    MPI.Printf0(stderr,"\t\tpressurize contours...\n");
    pressurize_contours();
    
    MPI.Printf0(stderr,"\t\tcompute gradP...\n");
    compute_gradP(MPI);
    
    MPI.Printf0(stderr,"\t\tcompute velocities...\n");
    compute_velocities();
}


#define __CHECK_CVG() do { if( Fabs(newP-oldP) > Fabs(ftol*newP) ) converged=0; } while(false)


int Workspace:: update_pressure(const mpi &MPI,
                                ColorType  c)
{
    
    //--------------------------------------------------------------------------
    // set the boundary conditions
    //--------------------------------------------------------------------------
    if(right_wall)
    {
        
    }
    else
    {
        for(unit_t j=outline.lower.y;j<=outline.upper.y;++j)
            P[j][upper.x] = P_user;
    }
    
    //--------------------------------------------------------------------------
    // prepare the pressure fields
    //--------------------------------------------------------------------------
    pressurize_bubbles();   //!< initialize Enter/Leave
    pressurize_contours();  //!< finalize   Enter/Leave according to bubbles
    
    //--------------------------------------------------------------------------
    // Evaluate gradient
    //--------------------------------------------------------------------------
    compute_gradP(MPI);
    
    //--------------------------------------------------------------------------
    // update according to color
    //--------------------------------------------------------------------------
    const unit_t shift[2] = { c, 1-c };
    assert(shift[0]==0||shift[0]==1);
    assert(shift[1]==0||shift[1]==1);
    int converged = 1;
    for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
    {
        const unit_t jm = j-1;
        const unit_t jp = j+1;
        for( unit_t i=bulk_imin + shift[j&1]; i <= bulk_imax; i +=2 )
        {
            if(B[j][i]<0)
            {
                const Real P_center  = P[j][i];
                const Real P_left    = L2[j][i-1].x;
                const Real P_right   = E2[j][i+1].x;
                const Real P_bottom  = L2[jm][i].y;
                const Real P_top     = E2[jp][i].y;
                const Real mid       = -(P_center+P_center);
                const Real Laplacian = (P_left+mid+P_right) * order2fac.x + (P_bottom+mid+P_top) * order2fac.y;
                const Real residue   = (Laplacian);
                const Real dP        = residue / W[j][i];
                
                P[j][i] -= dP;
                
                if( Fabs(dP) > Fabs(ftol*P[j][i]) )
                {
                    converged = 0;
                }
                
                DeltaP[j][i] = Laplacian;
                
            }
        }
    }
    
    //--------------------------------------------------------------------------
    // update the left boundary conditions
    //--------------------------------------------------------------------------
    const unit_t i0 = lower.x;
    const unit_t i1 = i0+1;
    const unit_t i2 = i1+1;
    
    for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
    {
        assert(B[j][i0]<0);
        if(B[j][i1]>=0)
        {
            // order 1 setting
            P[j][i0] = E1[j][i1].x;
        }
        else
        {
            const Real oldP = P[j][i0];
            const Real newP = P[j][i0] = (4.0 * P[j][i1] - E1[j][i2].x) / 3.0;
            __CHECK_CVG();
        }
    }
    
    //--------------------------------------------------------------------------
    // update the right boundary condition
    //--------------------------------------------------------------------------
    const unit_t in0 = upper.x;
    const unit_t in1 = in0-1;
    const unit_t in2 = in1-1;
    if( right_wall )
    {
        for(unit_t j=bulk_jmin;j<=bulk_jmax;++j)
        {
            assert(B[j][in0]<0);
            if(B[j][in1]>=0)
            {
                // order 1 setting
                P[j][in0] = L1[j][in1].x;
            }
            else
            {
                const Real oldP = P[j][in0];
                const Real newP = P[j][in0] = (4.0 * P[j][in1] - L1[j][in2].x) / 3.0;
                __CHECK_CVG();
            }
        }
    }
    
    
    //--------------------------------------------------------------------------
    // lower boundary condition
    //--------------------------------------------------------------------------
	if(MPI.IsFirst)
	{
		const unit_t j0 = lower.y;
		const unit_t j1 = j0+1;
		const unit_t j2 = j1+1;
        
		for(unit_t i=i1;i<upper.x;++i)
		{
			assert(B[j0][i]<0);
			if(B[j1][i]>=0)
			{
				// order 1 setting
				P[j0][i] = E1[j1][i].y;
			}
			else
			{
                const Real oldP = P[j0][i];
				const Real newP = P[j0][i] = (4.0 * P[j1][i] - E1[j2][i].y)/3.0;
                __CHECK_CVG();
			}
		}
        
        // lower left corner
        {
            const Real oldP = P[j0][i0];
            const Real newP = P[j0][i0] = (P[j1][i0] + P[j0][i1])/2;
            __CHECK_CVG();
        }
        
        if(right_wall)
        {
            //lower right corner
            const Real oldP = P[j0][in0];
            const Real newP = P[j0][in0] = (P[j1][in0]+P[j0][in1])/2;
            __CHECK_CVG();
        }
	}
	
    //--------------------------------------------------------------------------
    // upper boundary condition
    //--------------------------------------------------------------------------
	if(MPI.IsFinal)
	{
		const unit_t j0 = upper.y;
		const unit_t j1 = j0-1;
		const unit_t j2 = j1-1;
		for(unit_t i=i1;i<upper.x;++i)
		{
			assert(B[j0][i]<0);
			if(B[j1][i]>=0)
			{
				// order 1 setting
				P[j0][i] = L1[j1][i].y;
			}
			else
			{
                const Real oldP = P[j0][i];
				const Real newP = P[j0][i] = (4.0 * P[j1][i] - L1[j2][i].y)/3.0;
                __CHECK_CVG();
			}
		}
        
        {
            const Real oldP = P[j0][i0];
            const Real newP = P[j0][i0] = (P[j1][i0] + P[j0][i1])/2;
            __CHECK_CVG();
        }
        
        if( right_wall )
        {
            const Real oldP = P[j0][in0];
            const Real newP = P[j0][in0] = (P[j1][in0]+P[j0][in1])/2;
            __CHECK_CVG();
        }
	}
    
	sync1(MPI, P);
	
    return converged;
}
