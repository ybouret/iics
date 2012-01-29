#include "./domain2d.hpp"



using namespace Bubble;


static inline 
Real initRho( Real x, Real y )
{
	const Real r2 = x*x + y*y;
	if( sqrt(r2) <= 0.1 )
		return 1;
	else
		return 0.5; //+ 0.1 * ( 0.5 - alea<Real>() );
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
		mpi_next  = MPI.CommWorldNext();
		mpi_prev  = MPI.CommWorldPrev();
		MPI.Printf(stderr, "-- Rank %d | Size=%d | %d -> %d -> %d\n", mpi_rank, mpi_size, mpi_prev, mpi_rank, mpi_next );
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup  layouts
		//
		////////////////////////////////////////////////////////////////////////
		
		//======================================================================
		// geometry
		//======================================================================
		const Layout   full_layout( Coord(1,1), Coord(100,100) );
		vector<Layout> layouts(mpi_size,as_capacity);
		for( int r=0; r < mpi_size; ++r )
		{
			layouts.push_back( full_layout.split(mpi_rank,mpi_size) );
		}
		const Real   Lx=0.7,Ly=0.8;
		const Region full_region( Vertex(-Lx/2,-Ly/2), Vertex(Lx/2,Ly/2) );
		
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
		sim_ghosts.outer.lower.count = Coord(1,1);
		sim_ghosts.outer.lower.async = Coord(0,1);
		sim_ghosts.outer.lower.peers = Coord(-1,mpi_prev);
		
		sim_ghosts.outer.upper.count = Coord(1,1);
		sim_ghosts.outer.upper.async = Coord(0,1);
		sim_ghosts.outer.upper.peers = Coord(-1,mpi_next);
		
		sim_ghosts.inner.lower.count = Coord(1,1);
		sim_ghosts.inner.lower.async = Coord(0,1);
		sim_ghosts.inner.lower.peers = Coord(-1,mpi_prev);
		
		
		sim_ghosts.inner.upper.count = Coord(1,1);
		sim_ghosts.inner.upper.async = Coord(0,1);
		sim_ghosts.inner.upper.peers = Coord(-1,mpi_next);
		
		Domain domain( sim_layout, sim_ghosts, sim_region );
		MPI.Printf( stderr, "Rank %d: #asyncGhosts= %2lu, #plainGhosts= %2lu\n", mpi_rank, domain.async_ghosts, domain.plain_ghosts );
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// Initial Conditions
		//
		////////////////////////////////////////////////////////////////////////
		VdW fluid(1.1);
		{
			Fill::function2 f( cfunctor2(initRho) );
			Fill::with(f, domain["rho"], domain, domain.X, domain.Y );
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
			rwops<Real>::save_vtk( "fields2d_0.vtk", "", *pW, save_index, *pW );
		}
		
		
		for( int count=1; count <= 500; ++count )
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
				rwops<Real>::save_vtk( vformat("fields2d_%d.vtk",count), "", *pW, save_index, *pW );
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