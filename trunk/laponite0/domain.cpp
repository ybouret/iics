#include "./domain.hpp"

namespace Laponite
{
	
	const char * Domain::VarNames[] = 
	{
		"rho",
		"P",
		"Ux",
		"Uy",
		"dPdx",
		"dPdy"
	};
	
	const size_t Domain::VarCount = sizeof(VarNames) / sizeof(VarNames[0]);
	
	Domain:: ~Domain() throw()
	{
		
	}
	
	
	
	
	Domain:: Domain( const Layout &full_layout, const GhostsSetup &g, const Region &full_region ) :
	Workspace( full_layout.split(mpi_rank,mpi_size),
			  g,
			  Region::extract( full_region, full_layout, full_layout.split(mpi_rank,mpi_size) ),
			  VarCount,
			  VarNames),
	is_first( 0 == mpi_rank ),
	is_final( mpi_last == mpi_rank ),
	bulk( lower + Coord(0,is_first?1:0), upper-Coord(0,is_final?1:0)),
	requests( (2*async_ghosts) * 1 )
	{
		const mpi & MPI = mpi::instance();
		MPI.Printf( stderr, "Rank %d:  #requests= %lu\n", mpi_rank,  requests.count );
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
			for( size_t i = 1; i >0; --i )
			{
				
				MPI.Irecv( outer_ghost[i], outer_ghost.count, IICS_REAL, outer_ghost.peer, tag, MPI_COMM_WORLD, requests[ iRequest++ ] );
				
				inner_ghost.pull( Field[i], i );
				MPI.Isend( inner_ghost[i], inner_ghost.count, IICS_REAL, inner_ghost.peer, tag, MPI_COMM_WORLD, requests[ iRequest++ ] );
			}
			
		}
		
		//==============================================================
		// then copy the plain ghosts
		//==============================================================
		for( size_t g = plain_ghosts; g > 0; --g )
		{
			for( size_t i=1;i>0;--i)
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
	
	void Domain:: update_boundaries()
	{
		if( mpi_rank == 0 )
		{
			
		}
	}
	
	void Domain:: compute_P( const Layout &sub, const VanDerWaals &gas )
	{
		Workspace   &Field = *this;
		const Array &rho   = Field["rho"];
		Array       &P     = Field["P"];
		
		//-- first pass: compute P
		for( unit_t y=sub.upper.y; y >= sub.lower.y; --y )
		{
			for( unit_t x=sub.upper.x; x >= sub.lower.x; --x )
			{
				P[y][x] = gas.pressure( rho[y][x] );
			}
		}
		
#if 0
		//-- second pass: compute grapP and...
		Array &dPdx = Field["dPdx"];
		Array &dPdy = Field["dPdy"];
		
		for( unit_t y=sub.upper.y; y >= sub.lower.y; --y )
		{
			for( unit_t x=sub.upper.x; x >= sub.lower.x; --x )
			{
				dPdx[y][x] = (P[y][x+1] - P[y][x-1]);
				dPdy[y][x] = (P[y+1][x] - P[y-1][x]);
			}
		}
		
#endif
		
	}
	
	
}
