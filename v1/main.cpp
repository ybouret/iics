#include "common.hpp"

using namespace IICS;

static const char  *var_names[] = { "u", "v", "DeltaU","DeltaV" };
static const size_t var_count   = sizeof(var_names)/sizeof(var_names[0]);
static const size_t var_start   = 1; //-- variables from 1 to var_count


static inline 
void exchange( mpi &MPI, Workspace &W, const array<size_t> &variables_to_exchange )
{
	
	for( size_t g = W.ghosts; g>0; --g )
	{
		const Ghost &outer_ghost = W.outer_ghost(g);
		const Ghost &inner_ghost = W.inner_ghost(g);
		switch( inner_ghost.position )
		{
			case ghost_lower_x:
			case ghost_lower_y:
			case ghost_upper_x:
			case ghost_upper_y:
				for( size_t j=variables_to_exchange.size(); j > 0; --j )
				{
					const size_t iv = variables_to_exchange[j];
					Ghost::direct_copy( outer_ghost, inner_ghost, W[iv] );
				}
				break;
				
			default:
				break;
		}
	}
}


static inline Real FillSphere( Real x, Real y, Real z )
{
	if( x*x + y*y + z*z <= 1 )
	{
		return 1;
	}
	else {
		return 0;
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
		const size_t rank = MPI.CommWorldRank;
		const size_t size = MPI.CommWorldSize;
		
		std::cerr << "-- size = " << size << std::endl;
		std::cerr << "-- rank = " << rank << std::endl;
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup simulation: gather full layout and region
		//
		////////////////////////////////////////////////////////////////////////
		const Real Lx= 15.0, Ly=20.0, Lz=25.0;
		const Layout full_layout( Coord(1,1,1), Coord(10,15,20) );
		const Region full_region( Vertex(-Lx/2,-Ly/2,-Lz/2), Vertex(Lx/2,Ly/2,Lz/2) );
		MPI.Barrier(MPI_COMM_WORLD);
		if( rank == 0 )
		{
			std::cerr << "-- Full Layout: " << full_layout.lower << " -> " << full_layout.upper << " : width  = " << full_layout.width << std::endl;
			std::cerr << "-- Full Region: " << full_region.min   << " -> " << full_region.max   << " : length = " << full_region.length << std::endl;
		}
		
		////////////////////////////////////////////////////////////////////////
		//
		// split layout for this rank 
		//
		////////////////////////////////////////////////////////////////////////
		const Layout sim_layout = full_layout.split(rank, size);
		const Region sim_region( Region::extract( full_region, full_layout, sim_layout) );
		MPI.Barrier(MPI_COMM_WORLD);
		{
			std::cerr << "-- Rank " << rank << " Layout: " << sim_layout.lower << " -> " << sim_layout.upper << " : width  = " << sim_layout.width << std::endl;
			std::cerr << "-- Rank " << rank << " Region: " << sim_region.min   << " -> " << sim_region.max   << " : length = " << sim_region.length << std::endl;
		}
		
		////////////////////////////////////////////////////////////////////////
		//
		// create workspace with one ghost in each dimension for PBC and MPI
		//
		////////////////////////////////////////////////////////////////////////
		Workspace W(sim_layout.lower, sim_layout.upper, //-- indices
					Coord(1,1,1),     Coord(1,1,1),     //-- ghosts
					sim_region.min,   sim_region.max,   //-- space vertices
					var_start, var_count, var_names );  //-- variables
		const Layout &inside = W;
		vector<size_t> variables_to_exchange;
		variables_to_exchange.push_back( W("u") );
		variables_to_exchange.push_back( W("v") );
		
		////////////////////////////////////////////////////////////////////////
		//
		// Fill Workspace
		//
		////////////////////////////////////////////////////////////////////////
		FillFunctor fill_function( cfunctor3( FillSphere ) );
		Fill::with( fill_function, W["u"], inside, W.X, W.Y, W.Z ); //!< fill u layout (not the ghosts) with fill1
		W["v"].set_all( inside, 0.0 );
		
		////////////////////////////////////////////////////////////////////////
		//
		// Run
		//
		////////////////////////////////////////////////////////////////////////		
		Real dt = 0.001;
		
		for( size_t iter=1; iter <= 10; ++iter)
		{
			if( rank == 0 )
				std::cerr << "iter=" << iter << std::endl;
			exchange(MPI,W,variables_to_exchange);
			Laplacian::compute( W["DeltaU"], dt,  W["u"], W.inv_dsq, inside );
			Laplacian::compute( W["DeltaV"], dt,  W["v"], W.inv_dsq, inside );
			
			W["u"].add( W["DeltaU"], inside );
			W["v"].add( W["DeltaV"], inside );
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
