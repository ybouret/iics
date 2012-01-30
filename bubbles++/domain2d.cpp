#include "./domain2d.hpp"

namespace Bubble
{
	
	const char * Domain:: VarNames[] = { "rho", "U", "V", "P", "LU", "LV", "Lrho" };
	const size_t Domain:: VarCount   = sizeof(VarNames)/sizeof(VarNames[0]);
	const size_t Domain:: VarMPI     = 4;
	
	Domain:: ~Domain() throw() {}
	
	
	Domain:: Domain(const Layout      &L,
					const GhostsSetup &G,
					const Region      &R ) :
	Workspace(L,G,R,VarCount,VarNames),
	eta(0.001),
	gam(0.001),
	field_index(),
	chrono(),
	requests( VarMPI * (2*async_ghosts) )
	{
		const Workspace &self = *this;
		
		// register the variables to communicate
		for( size_t i=0; i < VarMPI; ++i )
		{
			field_index.push_back( self(VarNames[i]) );
		}
		
		// memory for async ghosts
		acquire_ghosts_data( field_index.size() );
		
	}
	
	
	double Domain:: exchanges_start()
	{
		chrono.start();
		_mpi::exchanges_start<Workspace>( *this, field_index, requests, 7 );
		
		return chrono.query();
	}
	
	double Domain:: exchanges_finish()
	{
		chrono.start();
		_mpi::exchanges_finish<Workspace>( *this, field_index, requests );
		
		return chrono.query();
	}
	
	static inline void pressure_cb( Real &P, const Real &rho, void *args )
	{
		const Fluid &fluid = *(Fluid *)args;
		P = fluid.pressure(rho);
	}
	
	double Domain:: compute_pressure( const Fluid &fluid )
	{
		chrono.start();
		Workspace &field = *this;
		void      *args = (void*)&fluid;
		field["P"].foreach( field["rho"], field, pressure_cb, args );
		return chrono.query();
	}
	
	static inline Real __laplacian( const Array &F, unit_t x, unit_t y, const Vertex &fac )
	{
		const Real F0 = F[y][x];
		const Real twoF0 = F0+F0;
		
		return fac.x * ( F[y][x+1] - twoF0 + F[y][x-1] ) + fac.y * ( F[y+1][x] - twoF0 + F[y-1][x] );
		
	}
	
	double Domain:: update_fields( Real dt )
	{
		chrono.start();
		Workspace &field = *this;
		Array &rho  = field["rho"];
		Array &U    = field["U"];
		Array &V    = field["V"];
		Array &P    = field["P"];
		Array &LU   = field["LU"];
		Array &LV   = field["LV"];
		Array &Lrho = field["Lrho"];
		
		const Real half = 0.5 * dt;
		const Real eta_fac = eta * dt;
		
		
		
		
		for( unit_t y=upper.y; y >= lower.y; --y )
		{
			for(unit_t x=upper.x; x >= lower.x; --x )
			{
				const Real divJ_x = inv_d.x * (U[y][x+1] - U[y][x-1]);
				const Real divJ_y = inv_d.y * (V[y+1][x] - V[y-1][x]);
				rho[y][x] -= half * ( divJ_x + divJ_y);
				
				LU[y][x]   = __laplacian(U,  x,y,inv_dsq);
				LV[y][x]   = __laplacian(V,  x,y,inv_dsq);
				Lrho[y][x] = __laplacian(rho,x,y,inv_dsq);
			}
		}
		
		for(unit_t x=upper.x+1; x >= lower.x-1; --x )
		{
			Lrho[upper.y+1][x] = __laplacian(rho,x,upper.y+1,inv_dsq);
			Lrho[lower.y-1][x] = __laplacian(rho,x,lower.y-1,inv_dsq);
		}
		
		for( unit_t y=upper.y; y >= lower.y; --y )
		{
			for(unit_t x=upper.x; x >= lower.x; --x )
			{
				U[y][x] += half * inv_d.x *( gam * (Lrho[y][x+1] - Lrho[y][x-1]) -  ( P[y][x+1] - P[y][x-1] ) ) + eta_fac * LU[y][x];
				V[y][x] += half * inv_d.y *( gam * (Lrho[y+1][x] - Lrho[y-1][x]) -  ( P[y+1][x] - P[y-1][x] ) ) + eta_fac * LV[y][x];
			}
		}
		
		return chrono.query();
	}
	
	
}
