#include "./domain.hpp"
#include "yocto/code/rand.hpp"

using namespace Laponite;

const Real        Lx = 1;
const Real        Ly = 2;
const Real        Lz = 0.01;

static inline
Real initRho( Real x, Real y )
{
	return 1e-3 * ( 0.1 + (1+sin(2*y)) + 0.05 * alea<Real>() );
}


static inline
double vproc( const Real &v )
{
	return v;
}

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
		mpi_last  = MPI.CommWorldLast;
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup domain info
		//
		////////////////////////////////////////////////////////////////////////
		
		//======================================================================
		// coordinates/space
		//======================================================================
		
		const Layout      full_layout( Coord(0,0),  Coord(99,199)   );
		const Region      full_region( Vertex(0,0), Vertex(Lx,Ly) );
		GhostsSetup       topology;
		
		//======================================================================
		//-- default topology
		//-- one ghost per dimension, async on Y, PBC on X
		//======================================================================
		topology.outer.lower.count = Coord(1,1);
		topology.outer.lower.async = Coord(0,1);
		topology.outer.lower.peers = Coord(-1,mpi_rank-1);
		
		topology.outer.upper.count = Coord(1,1);
		topology.outer.upper.async = Coord(0,1);
		topology.outer.upper.peers = Coord(-1,mpi_rank+1);
		
		topology.inner.lower.count = Coord(1,1);
		topology.inner.lower.async = Coord(0,1);
		topology.inner.lower.peers = Coord(-1,mpi_rank-1);
		
		topology.inner.upper.count = Coord(1,1);
		topology.inner.upper.async = Coord(0,1);
		topology.inner.upper.peers = Coord(-1,mpi_rank+1);
		
		
		//======================================================================
		// special case mpi_rank == 0
		//======================================================================
		if( mpi_rank == 0 )
		{
			//-- no lower outer ghost on Y
			topology.outer.lower.count = Coord(1,0);
			topology.outer.lower.async = Coord(0,0);
			topology.outer.lower.peers = Coord(-1,-1);
			
			//-- no lower inner ghost on Y
			topology.inner.lower.count = Coord(1,0);
			topology.inner.lower.async = Coord(0,0);
			topology.inner.lower.peers = Coord(-1,-1);
			
		}
		
		//======================================================================
		// special case mpi_rank == mpi_last
		//======================================================================
		if( mpi_rank == mpi_last )
		{
			//-- no upper outer ghost on Y
			topology.outer.upper.count = Coord(1,0);
			topology.outer.upper.async = Coord(0,0);
			topology.outer.upper.peers = Coord(-1,-1);
			
			//-- no upper inner ghost on Y
			topology.inner.upper.count = Coord(1,0);
			topology.inner.upper.async = Coord(0,0);
			topology.inner.upper.peers = Coord(-1,-1);
		}
		
		MPI.Printf( stderr, "Rank %d: %4ld -> %4d -> %4ld\n", mpi_rank, topology.outer.lower.peers.y, mpi_rank, topology.outer.upper.peers.y );
		MPI.Printf0(stderr, "\n");
		Domain domain(full_layout,topology,full_region);
		MPI.Printf( stderr, "Rank %d: #plainGhosts=%lu, #asyncGhosts=%lu\n", mpi_rank, domain.plain_ghosts, domain.async_ghosts);	
		
		
		//! argon, Lz thick at 25 C
		VanDerWaals gas( 1.363, 0.03219, 40e-3, Lz, 273.15 + 25 );
		
		{
			FillFunctor f( cfunctor2(initRho) );
			Fill::with( f, domain["rho"], domain, domain.X, domain.Y );
			//domain["rho"].set_all(domain,1);
		}
		
		domain.exchange_start(MPI);
		domain.exchange_finish(MPI);
		
		domain.compute_P(domain,gas);
		
		Real vmax, vmin;
		domain["rho"].get_min_max(vmin,NULL,vmax,NULL);
		domain["rho"].ppm(  vformat( "rho%d.ppm", mpi_rank ), "", domain, vproc, NULL, vmin, vmax );
		
		domain["P"].get_min_max(vmin,NULL,vmax,NULL);
		domain["P"].ppm(  vformat( "p%d.ppm", mpi_rank ), "", domain, vproc, NULL, vmin, vmax );
		
		
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
