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
	is_side( (0 == mpi_rank) || (mpi_last == mpi_rank) ),
	requests( int( is_side  ? 2 : 4 ) * number )
	{
		const mpi & MPI = mpi::instance();
		MPI.Printf( stderr, "Rank %d: is_side=%d, #requests= %lu\n", mpi_rank, is_side ?  1 : 0 , requests.count );
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
		const int  tag   = 100;
		//==============================================================
		// start async ghosts
		//==============================================================
		size_t iRequest = 0;
		for( size_t g = async_ghosts; g>0; --g )
		{
			const Ghost &outer_ghost = async_outer_ghost(g);
			const Ghost &inner_ghost = async_inner_ghost(g);
			for( size_t i = number; i >0; --i )
			{
				
				MPI.Irecv( outer_ghost[i], outer_ghost.count, IICS_REAL, outer_ghost.peer, tag, MPI_COMM_WORLD, requests[ iRequest++ ] );
			
				inner_ghost.pull( Field[i], i );
				MPI.Isend( inner_ghost[i], inner_ghost.count, IICS_REAL, inner_ghost.peer, tag, MPI_COMM_WORLD, requests[ iRequest++ ] );
			}
			
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
	
	
	void Domain:: exchange_finish( const mpi &MPI )
	{
		Workspace &Field = *this;
		MPI.Waitall( requests );
		for( size_t g = async_ghosts; g>0; --g )
		{
			const Ghost & inner_ghost = async_inner_ghost(g);
			for( size_t i=number;i>0;--i)
				inner_ghost.pull( Field[i], i );
		}
		
	}

	
}
