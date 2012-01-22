#include "common.hpp"
#include "yocto/cliff/rwops.hpp"
#include "yocto/auto-ptr.hpp"
#include "yocto/wtime.hpp"

using namespace IICS;

static const char  *var_names[] = { "u", "v", "Lu","Lv" };
static const size_t var_count   = sizeof(var_names)/sizeof(var_names[0]);
static const size_t var_start   = 1; //-- variables from 1 to var_count


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

static
inline void start_exchange( Workspace &W, const array<size_t> &variables, mpi::Requests &requests )
{
	static const mpi & MPI = mpi::instance();
	
	size_t iRequest = 0;
	
	for( size_t g = W.ghosts; g > 0; --g )
	{
		const Ghost &outer_ghost = W.outer_ghost(g); //!< outer ghost
		const Ghost &inner_ghost = W.inner_ghost(g); //!< corresponding inner ghost
		
		switch( outer_ghost.position )
		{
				
				//-- using ghosts for PBC
			case ghost_lower_x:
			case ghost_lower_y:
			case ghost_upper_x:
			case ghost_upper_y:
				assert( outer_ghost.is_async == false );
				assert( inner_ghost.is_async == false );
				for( size_t i=variables.size();i>0;--i)
				{
					Ghost::direct_copy( outer_ghost, inner_ghost, W[ variables[i] ] );
				}
				break;
				
				//-- using MPI
			case ghost_lower_z:
			case ghost_upper_z:
				assert( outer_ghost.is_async == true );
				assert( inner_ghost.is_async == true );
				assert( outer_ghost.peer >= 0 );
				assert( inner_ghost.peer >= 0 );
				
				for( size_t i= variables.size();i>0;--i)
				{
					MPI.Irecv( outer_ghost[i], outer_ghost.count, IICS_REAL, outer_ghost.peer, 0, MPI_COMM_WORLD, requests[ iRequest++] );
					
					inner_ghost.pull( W[ variables[i] ], i );
					MPI.Isend( inner_ghost[i], inner_ghost.count, IICS_REAL, inner_ghost.peer, 0, MPI_COMM_WORLD, requests[ iRequest++] );
				}
				break;
				
		}
		
		
	}
	
	
}


// compute laplacian on nucleus ( workspace layout minus the inner async ghosts )
static inline void start_laplacian( Workspace &W, Real dt, const array<size_t> &variables, const array<size_t> &laplacians )
{
	assert( variables.size() == laplacians.size() );
	const Layout &nucleus = W.nucleus;
	const Vertex &idsq    = W.inv_dsq;
	for( size_t i=variables.size();i>0;--i)
	{
		laplacian<Real,Real>::compute( W[ laplacians[i] ], dt, W[ variables[i] ], idsq, nucleus );
	}
}

// compute laplacian on inner async ghosts
static inline void finish_laplacian( Workspace &W, Real dt, const array<size_t> &variables, const array<size_t> &laplacians )
{
	assert( variables.size() == laplacians.size() );
	const Vertex &idsq    = W.inv_dsq;
	for( size_t g = W.async_ghosts; g>0;--g )
	{
		for( size_t i=variables.size();i>0;--i)
		{
			laplacian<Real,Real>::compute( W[ laplacians[i] ], dt, W[ variables[i] ], idsq, W.async_inner_ghost(g) );
		}
	}
}


static inline void finish_exchange(  Workspace &W, const array<size_t> &variables, mpi::Requests &requests )
{
	static const mpi & MPI = mpi::instance();
	MPI.Waitall(requests);
	for( size_t g = W.async_ghosts; g > 0; --g )
	{
		const Ghost &G = W.async_outer_ghost(g); 
		for( size_t i=variables.size();i>0;--i)
		{
			G.push( W[ variables[i] ], i );
		}
	}
}

//! add the dt x laplacian field to each associated variable
static inline void update_fields( Workspace &W,  const array<size_t> &variables, const array<size_t> &laplacians )
{
	assert( variables.size() == laplacians.size() );
	for( size_t i=variables.size();i>0;--i)
	{
		W[ variables[i] ].add( W[ laplacians[i] ], W );
	}
}

struct timings
{
	double total;
	double comm;
	double diff;
};

static inline 
void cycle(Workspace           &W, 
		   Real                 dt,
		   const array<size_t> &variables, 
		   const array<size_t> &laplacians, 
		   mpi::Requests       &requests,
		   wtime               &chrono,
		   timings             &tmx
		   )
{
	double t_comm   = 0;
	double t_diff   = 0;
	chrono.start();
	start_exchange(W, variables, requests );
	t_comm += chrono.query();
	
	chrono.start();
	start_laplacian(W, dt, variables, laplacians);
	t_diff += chrono.query();
	
	chrono.start();
	finish_exchange(W, variables, requests);
	t_comm += chrono.query();
	
	chrono.start();
	finish_laplacian(W, dt, variables, laplacians);
	update_fields(W, variables, laplacians);
	t_comm += chrono.query();
	
	tmx.total = (tmx.comm=t_comm) + (tmx.diff=t_diff);
}


int main( int argc, char *argv[] )
{
	try
	{
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup MPI
		//
		////////////////////////////////////////////////////////////////////////
		const mpi & MPI = mpi::init( &argc, &argv );
		rank  = MPI.CommWorldRank;
		size  = MPI.CommWorldSize;
		above = MPI.CommWorldNext();
		below = MPI.CommWorldPrev();
		MPI.Printf(stderr, "-- Rank %d | Size=%d | %d -> %d -> %d\n", rank, size, below, rank, above );
		
		
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup simulation: create full layout and region
		//
		////////////////////////////////////////////////////////////////////////
		const Real Lx= 15.0, Ly=20.0, Lz=25.0;
		const Layout full_layout( Coord(1,1,1), Coord(10,15,20) );
		const Region full_region( Vertex(-Lx/2,-Ly/2,-Lz/2), Vertex(Lx/2,Ly/2,Lz/2) );
		
		MPI.Printf0( stderr, "-- Full Layout: [%ld %ld %ld] -> [%ld %ld %ld] :  width=[ %ld %ld %ld ]\n", 
					full_layout.lower.x,
					full_layout.lower.y,
					full_layout.lower.z,
					full_layout.upper.x,
					full_layout.upper.y,
					full_layout.upper.z,
					full_layout.width.x,
					full_layout.width.y,
					full_layout.width.z
					);
		MPI.Printf0( stderr, "-- Full Region: : [%.2g %.2g %.2g] -> [%.2g %.2g %.2g] :  length=[ %.2g %.2g %.2g ]\n\n",
					full_region.min.x,
					full_region.min.y,
					full_region.min.z,
					full_region.max.x,
					full_region.max.y,
					full_region.max.z,
					full_region.length.x,
					full_region.length.y,
					full_region.length.z );
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup simulation layout and region: split along last dimesion
		//
		////////////////////////////////////////////////////////////////////////
		const Layout sim_layout( full_layout.split(rank, size) );
		MPI.Printf( stderr, "-- Rank %d Layout: [%ld %ld %ld] -> [%ld %ld %ld] :  width=[ %ld %ld %ld ]\n", 
				   rank,
				   sim_layout.lower.x,
				   sim_layout.lower.y,
				   sim_layout.lower.z,
				   sim_layout.upper.x,
				   sim_layout.upper.y,
				   sim_layout.upper.z,
				   sim_layout.width.x,
				   sim_layout.width.y,
				   sim_layout.width.z
				   );
		const Region sim_region( Region::extract( full_region, full_layout, sim_layout ) );
		MPI.Printf( stderr, "-- Rank %d Region: : [%.2g %.2g %.2g] -> [%.2g %.2g %.2g] :  length=[ %.2g %.2g %.2g ]\n",
				   rank,
				   sim_region.min.x,
				   sim_region.min.y,
				   sim_region.min.z,
				   sim_region.max.x,
				   sim_region.max.y,
				   sim_region.max.z,
				   sim_region.length.x,
				   sim_region.length.y,
				   sim_region.length.z );
		
		////////////////////////////////////////////////////////////////////////
		//
		// create ghosts setup: 1 in each dimension for PBC
		// the ghost in the last dimension is deferred for I/O overlapping
		//
		////////////////////////////////////////////////////////////////////////
		GhostsInfos ghosts_up_and_lo( Coord(1,1,1), Coord(0,0,1) );
		GhostsSetup sim_ghosts(ghosts_up_and_lo,ghosts_up_and_lo);
		
		////////////////////////////////////////////////////////////////////////
		//
		// create the workspace and assign async ghosts peers
		//
		////////////////////////////////////////////////////////////////////////
		Workspace W( sim_layout, sim_ghosts, sim_region, var_start, var_count, var_names );
		setup_ghosts_peers( W );
		
		const Layout &inside  = W;           //!< inside layout, without ghosts
		MPI.Printf0( stderr, "\n");
		
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// prepare the variables
		//
		////////////////////////////////////////////////////////////////////////
		vector<size_t> variables;
		vector<size_t> laplacians;
		variables.push_back( W("u") ); laplacians.push_back( W("Lu") );
		variables.push_back( W("v") ); laplacians.push_back( W("Lv") );
		
		const size_t num_fields = variables.size();
		
		////////////////////////////////////////////////////////////////////////
		//
		// create MPI requests and ghosts memory
		//
		////////////////////////////////////////////////////////////////////////
		mpi::Requests requests( 4*num_fields );
		W.acquire_ghosts_data( num_fields );
		
		
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// Initialize fields
		//
		////////////////////////////////////////////////////////////////////////
		for( size_t i=num_fields; i>0; --i )
		{
			W[ variables[i] ].set_all( inside, 1+rank );
		}
		
		
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// Run...
		//
		////////////////////////////////////////////////////////////////////////
		wtime chrono;
		const Real dt = 1e-4;
		timings    tmx;
		
		for( size_t iter=1; iter <= 100; ++iter )
		{
			cycle(W, dt, variables, laplacians, requests, chrono, tmx);
			MPI.Printf0( stderr, "time: %8.3f comm=%8.3f diff=%8.3f\n", 1000 * tmx.total, 1000 * tmx.comm, 1000 * tmx.diff);
		}
		
		
		return 0;
	}
	catch( const exception &e )
	{
		std::cerr << std::endl;
		std::cerr << "******** " << e.what() << std::endl;
		std::cerr << "******** " << e.when() << std::endl;
	}
	catch(...)
	{
		std::cerr << std::endl;
		std::cerr << "******** Unhandled Exception!" << std::endl;
	}
	
	return -1;
}
