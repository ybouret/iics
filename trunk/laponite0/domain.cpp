#include "./domain.hpp"

namespace Laponite
{
	
	const char * Domain::VarNames[] = 
	{
		"rho",
		"P",
		"Ux",
		"Uy"
	};
	
	const size_t Domain::VarCount = sizeof(VarNames) / sizeof(VarNames[0]);
	
	Domain:: ~Domain() throw()
	{
		
	}
	
	
	static inline void __get_peer_for( const Ghost &G )
	{
		
		assert( G.is_async );
		switch( G.position )
		{
			case ghost_lower_y:
				assert(mpi_prev>=0);
				G.peer = mpi_prev;
				return;
				
			case ghost_upper_y:
				assert(mpi_next>=0);
				G.peer = mpi_next;
				return;
				
			default:
				throw exception( "No peer for for ghost@%s", G.label() );
		}
	}
	
	
	Domain:: Domain( const Layout &full_layout, const GhostsSetup &g, const Region &full_region ) :
	Workspace( full_layout.split(mpi_rank,mpi_size),
			  g,
			  Region::extract( full_region, full_layout, full_layout.split(mpi_rank,mpi_size) ),
			  VarCount,
			  VarNames),
	requests( number )
	{
		
		//======================================================================
		// set peer
		//======================================================================
		for( size_t g = async_ghosts; g >0; --g )
		{
			__get_peer_for( async_outer_ghost(g) );
			__get_peer_for( async_inner_ghost(g) );
		}
		
		//======================================================================
		// data for async ghosts
		//======================================================================
		acquire_ghosts_data( number );
		
		
		
	}
	
	void Domain:: exchange_start( const mpi &MPI )
	{
		Workspace &Field = *this;
		//==============================================================
		// start async ghosts
		//==============================================================
		for( size_t g = async_ghosts; g>0; --g )
		{
			const Ghost &outer_ghost = async_outer_ghost(g);
		}
		
		//==============================================================
		// then plain ghosts
		//==============================================================
		for( size_t g = plain_ghosts; g > 0; --g )
		{
			for( size_t i=number;i>0;--i)
			{
				Ghost::direct_copy( plain_outer_ghost(g), plain_inner_ghost(g), Field[i] );
			}
		}
		
	}
	
	
}
