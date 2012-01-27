#include "./domain.hpp"

namespace Bubble
{
	
	const char * Domain:: VarNames[] = { "rho", "U", "V" };
	const size_t Domain:: VarCount   = sizeof(VarNames)/sizeof(VarNames[0]);
	
	static const size_t __mpi_count  = 3;
	
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
		_mpi<Workspace>::exchanges_start( *this, field_index, requests, 7 );
		
		return chrono.query();
	}
	
	double Domain:: exchanges_finish()
	{
		chrono.start();
		_mpi<Workspace>::exchanges_finish( *this, field_index, requests );
		
		return chrono.query();
	}
	
}
