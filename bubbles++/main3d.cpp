#include "./domain3d.hpp"



using namespace Bubble;


static inline 
Real initRho( Real x, Real y, Real z )
{
	const Real r2 = x*x + y*y + z*z;
	if( sqrt(r2) <= 0.1 )
		return 1;
	else
		return 0.5;// + 0.1 * ( 0.5 - alea<Real>() );
}

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
		const Layout   full_layout( Coord(1,1,1), Coord(20,30,40) );
		vector<Layout> layouts(mpi_size,as_capacity);
		for( int r=0; r < mpi_size; ++r )
		{
			layouts.push_back( full_layout.split(mpi_rank,mpi_size) );
		}
		const Real Lx=0.7,Ly=0.8,Lz=0.9;
		const Region full_region( Vertex(-Lx/2,-Ly/2,-Lz/2), Vertex(Lx/2,Ly/2,Lz/2) );
		
		const Layout &sim_layout = layouts[ mpi_rank + 1 ];
		const Region &sim_region = Region::extract( full_region, full_layout, sim_layout );
		GhostsSetup   sim_ghosts;
		
		//======================================================================
		// one full workspace
		//======================================================================
		auto_ptr<Workspace> pW;
		vector<size_t>      save_index;
		Array              *full_rho = NULL;
		Array              *full_P   = NULL;
		if( 0 == mpi_rank )
		{
			const char *fullNames[] = { "rho", "P" };
			pW.reset( new Workspace(full_layout,sim_ghosts,full_region,sizeof(fullNames)/sizeof(fullNames[0]), fullNames) );
			save_index.push_back(1);
			save_index.push_back(2);
			full_rho = & (*pW)[ "rho" ];
			full_P   = & (*pW)[ "P"   ];
		}
		
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
		VdW fluid(1.1);
		{
			Fill::function3 f( cfunctor3(initRho) );
			Fill::with(f, domain["rho"], domain, domain.X, domain.Y, domain.Z );
			domain.compute_pressure(fluid);
		}
		
		////////////////////////////////////////////////////////////////////////
		//
		// simulation
		//
		////////////////////////////////////////////////////////////////////////
		
		Real dt = 0.001;
		_mpi::collect0<Real>( full_rho, domain["rho"], full_layout );
		_mpi::collect0<Real>( full_P,   domain["P"],   full_layout );
		
		if( 0 == mpi_rank )
		{
			rwops<Real>::save_vtk( "fields3d_0.vtk", "", *pW, save_index, *pW );
		}
		
		
		for( int count=1; count <= 100; ++count )
		{
			MPI.Printf0( stderr, "count=%5d      \r", count );
			for( size_t iter=0; iter <10; ++iter )
			{
				domain.exchanges_start();
				domain.exchanges_finish();
				domain.update_fields(dt);
				domain.compute_pressure(fluid);
			}
			
			_mpi::collect0<Real>( full_rho, domain["rho"], full_layout );
			_mpi::collect0<Real>( full_P,   domain["P"],   full_layout );
			if( 0 == mpi_rank )
			{
				rwops<Real>::save_vtk( vformat("fields3d_%d.vtk",count), "", *pW, save_index, *pW );
			}
			
		}
		MPI.Printf0( stderr, "\n");
		
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