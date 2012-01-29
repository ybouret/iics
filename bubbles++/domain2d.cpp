#include "./domain2d.hpp"

namespace Bubble
{
	
	const char * Domain:: VarNames[] = { "rho", "U", "V", "P" };
	const size_t Domain:: VarCount   = sizeof(VarNames)/sizeof(VarNames[0]);
	const size_t Domain:: VarMPI     = 4;
	
	Domain:: ~Domain() throw() {}
	
	
	Domain:: Domain(const Layout      &L,
					const GhostsSetup &G,
					const Region      &R ) :
	Workspace(L,G,R,VarCount,VarNames),
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
	
	double Domain:: update_fields( Real dt )
	{
		chrono.start();
		Workspace &field = *this;
		Array &rho = field["rho"];
		Array &U   = field["U"];
		Array &V   = field["V"];
		Array &P   = field["P"];
		const Real half = 0.5 * dt;
		
		for( unit_t y=upper.y; y >= lower.y; --y )
		{
			for(unit_t x=upper.x; x >= lower.x; --x )
			{
				const Real divJ_x = inv_d.x * (U[y][x+1] - U[y][x-1]);
				const Real divJ_y = inv_d.y * (V[y+1][x] - V[y-1][x]);
				rho[y][x] -= half * ( divJ_x + divJ_y);
			}
		}
		
		for( unit_t y=upper.y; y >= lower.y; --y )
		{
			for(unit_t x=upper.x; x >= lower.x; --x )
			{
				U[y][x] -= half * ( P[y][x+1] - P[y][x-1] );
				V[y][x] -= half * ( P[y+1][x] - P[y-1][x] );
			}
		}
		
		return chrono.query();
	}
	
	
}
