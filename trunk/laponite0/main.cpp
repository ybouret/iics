#include "./domain.hpp"
#include "yocto/code/rand.hpp"

using namespace Laponite;

int main( int argc, char *argv[] )
{
	
	try
	{
		_rand.wseed();
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup MPI and ghosts
		//
		////////////////////////////////////////////////////////////////////////
		const mpi & MPI = mpi::init( &argc, &argv );
		mpi_rank  = MPI.CommWorldRank;
		mpi_size  = MPI.CommWorldSize;
	
		GhostsInfos  ghosts_lo;
		GhostsInfos  ghosts_up;
		
		
		if( mpi_size > 1 )
		{
				if( mpi_rank > 0 )
				{
					mpi_prev = mpi_rank -1;
					ghosts_lo.count = Coord(1,1); //! one lower ghost in X and Y
					ghosts_lo.async = Coord(0,1); //! communicating in Y direction
				}
				else
				{
					mpi_prev = MPI.CommWorldRankMax;
					ghosts_lo.count = Coord(1,0); //!< only on X for rank=0
					ghosts_lo.async = Coord(0,0); //!< plain ghosts for X
				}
			
				if( mpi_rank < MPI.CommWorldRankMax )
				{
					mpi_next = mpi_rank + 1;
					ghosts_up.count = Coord(1,1); //! one upper ghost in X and Y
					ghosts_up.async = Coord(0,1); //! communicating in Y direction
				}
				else
				{	
					mpi_next = 0;
					ghosts_up.count = Coord(1,0); //! only on X for rank=rank_max
					ghosts_up.async = Coord(0,0); //! plain ghosts for X
				}
		}
		else
		{
			ghosts_lo.count = ghosts_up.count = Coord(1,0); //! only on X
			ghosts_lo.async = ghosts_up.async = Coord(0,0); //! only 
		}
		
		MPI.Printf( stderr, "Rank %d: %d -> %d -> %d | ghosts_lo: (%ld,%ld) | ghosts_up (%ld,%ld)\n",
				   mpi_rank, mpi_prev, mpi_rank, mpi_next,
				   ghosts_lo.count.x, ghosts_lo.count.y,
				   ghosts_up.count.x, ghosts_up.count.y);
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup domain info
		//
		////////////////////////////////////////////////////////////////////////
		const Real        Lx = 1;
		const Real        Ly = 2;
		const Layout      full_layout( Coord(0,0),  Coord(9,14)   );
		const Region      full_region( Vertex(0,0), Vertex(Lx,Ly) );
		const GhostsSetup topology( ghosts_lo, ghosts_up );
		
		Domain domain(full_layout,topology,full_region);
		MPI.Printf( stderr, "Rank %d: #plainGhosts=%lu, #asyncGhosts=%lu\n", mpi_rank, domain.plain_ghosts, domain.async_ghosts);	
		
		
			
		return 0;
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