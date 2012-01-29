#include "./domain.hpp"

namespace Bubble
{
	
	const char * Domain:: VarNames[] = { "rho", "U", "V", "W", "P" };
	const size_t Domain:: VarCount   = sizeof(VarNames)/sizeof(VarNames[0]);
	
	static const size_t __mpi_count  = 5;
	
	Domain:: ~Domain() throw() {}
	
	
	Domain:: Domain(const Layout      &L,
					const GhostsSetup &G,
					const Region      &R ) :
	Workspace(L,G,R,VarCount,VarNames),
	field_index(),
	chrono(),
	requests( __mpi_count * (2*async_ghosts) )
	{
		const Workspace &self = *this;
		
		// register the variables to communicate
		for( size_t i=0; i < __mpi_count; ++i )
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
		Array &W   = field["W"];
		Array &P   = field["P"];
		const Real half = 0.5 * dt;
		for( unit_t z=upper.z; z >= lower.z; --z )
		{
			for( unit_t y=upper.y; y >= lower.y; --y )
			{
				for(unit_t x=upper.x; x >= lower.x; --x )
				{
					const Real divJ_x = inv_d.x * (U[z][y][x+1] - U[z][y][x-1]);
					const Real divJ_y = inv_d.y * (V[z][y+1][x] - V[z][y-1][x]);
					const Real divJ_z = inv_d.z * (W[z+1][y][x] - W[z-1][y][x]);
					rho[z][y][x] -= half * ( divJ_x + divJ_y + divJ_z);
				}
			}
		}
		
		for( unit_t z=upper.z; z >= lower.z; --z )
		{
			for( unit_t y=upper.y; y >= lower.y; --y )
			{
				for(unit_t x=upper.x; x >= lower.x; --x )
				{
					U[z][y][x] -= half * ( P[z][y][x+1] - P[z][y][x-1] );
					V[z][y][x] -= half * ( P[z][y+1][x] - P[z][y-1][x] );
					W[z][y][x] -= half * ( P[z+1][y][x] - P[z-1][y][x] );
				}
			}
		}
		
		return chrono.query();
	}

	
}
