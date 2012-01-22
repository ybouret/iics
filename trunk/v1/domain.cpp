#include "domain.hpp"


namespace IICS
{
	
	Domain:: ~Domain() throw()
	{
	}
	
	
	
	
	const Coord       Domain::GCount(1,1,1);
	const Coord       Domain::GAsync(0,0,1);
	const GhostsInfos Domain::GInfos( Domain::GCount, Domain::GAsync );
	const GhostsSetup Domain::GSetup( GInfos, GInfos );
	
	
	static inline void __get_peer_for( const Ghost &G )
	{
		assert(above>=0);
		assert(below>=0);
		switch( G.position )
		{
			case ghost_lower_z:
				G.peer =  below;
				return;
				
			case ghost_upper_z:
				G.peer = above;
				return;
				
			default:
				throw exception( "No peer for for ghost@%s", G.label() );
		}
	}
	
	static inline void setup_ghosts_peers( Workspace &W )
	{
		for( size_t g = W.async_ghosts; g>0; --g )
		{
			__get_peer_for( W.async_inner_ghost(g) );
			__get_peer_for( W.async_outer_ghost(g) );
		}
	}
	
	
	
	Domain:: Domain( const Layout &full_layout,
					const Region &full_region,
					size_t        fields,
					const char   *names[]) :
	Workspace(full_layout.split( rank, size ),
			  GSetup,
			  Region::extract( full_region, full_layout, full_layout.split( rank, size ) ),
			  VarStart,
			  fields*2,
			  names ),
	field_index(fields,as_capacity),
	delta_index(fields,as_capacity),
	num_fields(fields),
	requests( 4 * num_fields )
	{
		
		//! compute variables and laplacians indices
		for( size_t i=0; i < fields; ++i )
		{
			const size_t j = VarStart+i;
			field_index.push_back( j );
			delta_index.push_back( j+fields );
		}
		
		//! prepare ghosts data
		acquire_ghosts_data( num_fields );
		
		//! prepare ghosts connectivity
		setup_ghosts_peers(*this);
		
		
	}
	
}
