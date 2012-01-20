#include "common.hpp"

using namespace IICS;

static const char  *var_names[] = { "u", "v", "DeltaU","DeltaV" };
static const size_t var_count   = sizeof(var_names)/sizeof(var_names[0]);
static const size_t var_start   = 1; //-- variables from 1 to var_count


static
inline void start_exchange( Workspace &W, const array<size_t> &var )
{
	for( size_t g = W.ghosts; g > 0; --g )
	{
		
	}
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
		mpi & MPI = mpi::init( &argc, &argv );
		const int rank = MPI.CommWorldRank;
		const int size = MPI.CommWorldSize;
		
		MPI.Printf(stderr, "-- Rank %d | Size=%d\n", rank, size );
		
		
		
		
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
		GhostsInfo  ghosts_up_and_lo( Coord(1,1,1), Coord(0,0,1) );
		GhostsSetup sim_ghosts(ghosts_up_and_lo,ghosts_up_and_lo);
		
		////////////////////////////////////////////////////////////////////////
		//
		// create the workspace
		//
		////////////////////////////////////////////////////////////////////////
		Workspace W( sim_layout, sim_ghosts, sim_region, var_start, var_count, var_names );
		//const Layout &inside  = W;           //!< inside layout, without ghosts
		//const Layout &nucleus = W.nucleus;   //!< 
		MPI.Printf0( stderr, "\n");
		MPI.Printf( stderr, "Rank %d: #ghosts= %lu\n", rank, W.ghosts );
		MPI.Printf( stderr, "Rank %d: #deferred_ghost=%lu\n", rank, W.deferred_ghosts);
		assert( W.deferred_ghosts == 2 );
		
		
		
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// prepare the variables
		//
		////////////////////////////////////////////////////////////////////////
		vector<size_t> cid; //!< component ID
		cid.push_back( W("u") );
		cid.push_back( W("v") );
		
		const size_t num_cid = cid.size();
		
		////////////////////////////////////////////////////////////////////////
		//
		// create MPI requests
		//
		////////////////////////////////////////////////////////////////////////
		const size_t num_requests = 4 * num_cid;
		
		
		start_exchange( W, cid );
		
		
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
