#include "./domain.hpp"

#include "yocto/code/rand.hpp"

using namespace Bubble;

int main( int argc, char *argv[] )
{
	
	try
	{
		_rand.wseed();
		////////////////////////////////////////////////////////////////////////
		//
		// setup MPI
		//
		////////////////////////////////////////////////////////////////////////
		const mpi & MPI = mpi::init( &argc, &argv );
		mpi_rank  = MPI.CommWorldRank;
		mpi_size  = MPI.CommWorldSize;
		mpi_last  = MPI.CommWorldLast;
		mpi_above = MPI.CommWorldNext();
		mpi_below = MPI.CommWorldPrev();
		MPI.Printf(stderr, "-- Rank %d | Size=%d | %d -> %d -> %d\n", mpi_rank, mpi_size, mpi_below, mpi_rank, mpi_above );
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup  layouts
		//
		////////////////////////////////////////////////////////////////////////
		
		//======================================================================
		// geometry
		//======================================================================
		const Layout   full_layout( Coord(1,1,1), Coord(10,15,20) );
		vector<Layout> layouts(mpi_size,as_capacity);
		for( int r=0; r < mpi_size; ++r )
		{
			layouts.push_back( full_layout.split(mpi_rank,mpi_size) );
		}
		const Region full_region( Vertex(0,0,0), Vertex(0.7,0.8,0.9) );
		
		const Layout &sim_layout = layouts[ mpi_rank + 1 ];
		const Region &sim_region = Region::extract( full_region, full_layout, sim_layout );
		GhostsSetup   sim_ghosts;
		
		//======================================================================
		// topology: 1 ghost in each direction: x,y : PBC, z: MPI
		//======================================================================
		sim_ghosts.outer.lower.count = Coord(1,1,1);
		sim_ghosts.outer.lower.async = Coord(0,0,1);
		sim_ghosts.outer.lower.peers = Coord(-1,-1,mpi_below);
		
		sim_ghosts.outer.upper.count = Coord(1,1,1);
		sim_ghosts.outer.upper.async = Coord(0,0,1);
		sim_ghosts.outer.upper.peers = Coord(-1,-1,mpi_above);
		
		sim_ghosts.inner.lower.count = Coord(1,1,1);
		sim_ghosts.inner.lower.async = Coord(0,0,1);
		sim_ghosts.inner.lower.peers = Coord(-1,-1,mpi_below);

		
		sim_ghosts.inner.upper.count = Coord(1,1,1);
		sim_ghosts.inner.upper.async = Coord(0,0,1);
		sim_ghosts.inner.upper.peers = Coord(-1,-1,mpi_above);
		
		Domain domain( sim_layout, sim_ghosts, sim_region );
		MPI.Printf( stderr, "Rank %d: #asyncGhosts= %2lu, #plainGhosts= %2lu\n", mpi_rank, domain.async_ghosts, domain.plain_ghosts );
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// Initial Conditions
		//
		////////////////////////////////////////////////////////////////////////
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// simulation
		//
		////////////////////////////////////////////////////////////////////////
		domain.exchanges_start();
		domain.exchanges_finish();
		
	}
	catch( const exception &e )
	{
		std::cerr.flush();
		std::cerr << std::endl;
		std::cerr << "******** " << e.what() << std::endl;
		std::cerr << "******** " << e.when() << std::endl;
	}
	catch(...)
	{
		std::cerr.flush();
		std::cerr << std::endl;
		std::cerr << "******** Unhandled Exception!" << std::endl;
	}
	
	return -1;
	
}