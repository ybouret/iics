#include "domain.hpp"


namespace IICS
{
	
	Domain:: ~Domain() throw()
	{
	}
	
	
	
	
	Domain:: Domain(const Layout      &full_layout,
					const GhostsSetup &setup,
					const Region      &full_region,
					size_t             fields,
					const char        *names[]) :
	Workspace(full_layout.split( mpi_rank, mpi_size ),
			  setup,
			  Region::extract( full_region, full_layout, full_layout.split(mpi_rank,mpi_size) ),
			  fields*2+1,
			  names ),
	field_index(fields,as_capacity),
	delta_index(fields,as_capacity),
	num_fields(fields),
	requests( 4 * num_fields ),
	chrono(),
	odeint(),
	var( num_fields, 0 ),
	irx()
	{
		
		//! compute variables and laplacians indices
		for( size_t i=1; i <= fields; ++i )
		{
			field_index.push_back( i );
			delta_index.push_back( i+fields );
		}
		
		//! prepare ghosts data
		acquire_ghosts_data( num_fields );
		
		//! prepare ghosts connectivity
		//setup_ghosts_peers(*this);
		
		//! prepare driver
		odeint.start( num_fields );
		
		//! offsets to perform reaction
		Workspace &Field = *this;
		Field[1].load_offsets( irx, Field[1], outline );
		
	}
	
	double Domain:: start_exchanges( const mpi &MPI )
	{
		chrono.start();
		Workspace &Field  = *this;
		
		//-- start asynchronous exchanges
		size_t iRequest = 0;
		for( size_t g = async_ghosts; g>0; --g )
		{
			const Ghost &outer_ghost = async_outer_ghost(g);
			const Ghost &inner_ghost = async_inner_ghost(g);
			for( size_t i= num_fields;i>0;--i)
			{
				MPI.Irecv( outer_ghost[i], outer_ghost.count, IICS_REAL, outer_ghost.peer, 0, MPI_COMM_WORLD, requests[ iRequest++] );
				
				inner_ghost.pull( Field[ field_index[i] ], i );
				MPI.Isend( inner_ghost[i], inner_ghost.count, IICS_REAL, inner_ghost.peer, 0, MPI_COMM_WORLD, requests[ iRequest++] );
			}
		}
		
		//-- perform plain exchanges
		for( size_t g = plain_ghosts; g>0; --g )
		{
			const Ghost &outer_ghost = plain_outer_ghost(g);
			const Ghost &inner_ghost = plain_inner_ghost(g);
			for( size_t i= num_fields;i>0;--i)
			{
				Ghost::direct_copy( outer_ghost, inner_ghost, Field[ field_index[i] ] );
			}
			
		}
		
		return chrono.query();
		
	}
	
	double Domain:: finish_exchanges(  const mpi &MPI )	
	{
		chrono.start();
		Workspace &Field  = *this;
		MPI.Waitall(requests);
		for( size_t g = async_ghosts; g > 0; --g )
		{
			const Ghost &G = async_outer_ghost(g); 
			for( size_t i= num_fields;i>0;--i)
			{
				G.push( Field[ field_index[i] ], i );
			}
		}
		
		return chrono.query();
	}
	
	double Domain:: start_laplacian( Real dt )
	{
		chrono.start();
		Workspace &Field  = *this;
		for( size_t i=num_fields;i>0;--i)
		{
			laplacian<Real,Real>::compute( Field[ delta_index[i] ], dt, Field[ field_index[i] ], inv_dsq, nucleus );
		}
		return chrono.query();
	}
	
	double Domain:: finish_laplacian( Real dt )
	{
		chrono.start();
		Workspace &Field  = *this;
		for( size_t g = async_ghosts; g>0;--g )
		{
			for( size_t i=num_fields;i>0;--i)
			{
				laplacian<Real,Real>::compute( Field[ delta_index[i] ], dt, Field[ field_index[i] ], inv_dsq, async_inner_ghost(g) );
			}
		}
		return chrono.query();
	}
	
	double Domain:: update()
	{
		chrono.start();
		Workspace &Field  = *this;
		for( size_t i= num_fields;i>0;--i)
		{
			Field[ field_index[i] ].add( Field[ delta_index[i] ], *this );
		}
		return chrono.query();
	}
	
	double Domain:: reaction( ODE_Function &F, double t_curr, double t_next )
	{
		chrono.start();
		Workspace &Field  = *this;
		Real      *h      = Field["h"].entry;
		
		for( size_t i = irx.size(); i >0; --i )
		{
			//-- get the offset
			const size_t j = irx[i];
			
			//-- query data
			query( var, field_index, j);
			
			//-- solve
			odeint( F, var, t_curr, t_next, h[j] );
			
			//-- store data back
			store( var, field_index, j );
			
		}
		return chrono.query();
	}
	
	
	void Domain::cycle( Real dt, const mpi &MPI, Timings &timings )
	{
		timings.t_comm += start_exchanges(MPI);
		timings.t_diff += start_laplacian(dt);
		timings.t_comm += finish_exchanges(MPI);
		timings.t_diff += finish_laplacian(dt);
		timings.t_diff += update();
	}
	
	
}
