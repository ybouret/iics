#include "domain.hpp"
#include "reaction.hpp"

#include "yocto/cliff/rwops.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/ios/icstream.hpp"

#include "yocto/auto-ptr.hpp"
#include "yocto/filesys/local-fs.hpp"
#include "yocto/string/vfs-utils.hpp"

#include "yocto/code/rand.hpp"

using namespace IICS;
using namespace filesys;


////////////////////////////////////////////////////////////////////////////////
//
// I/O routines to collect sub domains
//
////////////////////////////////////////////////////////////////////////////////
static const char resdir[] = "results";

static inline void save_all( size_t iter, auto_ptr<Workspace> &pW0, const Workspace &W )
{
	static const mpi & MPI = mpi::instance();
	vfs &fs = local_fs::instance();
	const string outdir = _vfs::to_directory( resdir );
	fs.create_dir( outdir, true);
	
	{
		const string   output_name( outdir + vformat("sub%d.dat", mpi_rank) );
		ios::ocstream  output( output_name, false );
		W[1].save( output, W);
	}
	
	MPI.Barrier( MPI_COMM_WORLD );
	if( 0 == mpi_rank )
	{
		Workspace &W0 = *pW0;
		Array     &A  = W0[1];
		for( int r=0; r < mpi_size; ++r )
		{
			const string  input_name( outdir + vformat("sub%d.dat",r) );
			ios::icstream input( input_name );
			const Layout  sub( W0.split(r,mpi_size) );
			A.load( input, sub );
		}
		const string  var_name = W.name(1);
		const string  output_name( outdir + var_name + vformat("%lu.vtk",iter) );
		rwops<Real>::save_vtk( output_name, "", var_name, A, A, W0.region.min, W0.delta);		
	}
	MPI.Barrier( MPI_COMM_WORLD );
}

////////////////////////////////////////////////////////////////////////////////
//
// Initializing functions
//
////////////////////////////////////////////////////////////////////////////////
static inline
Real initV( Real x, Real y, Real z )
{
	
	return 1 + (0.5 - alea<Real>());
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
		mpi_above = MPI.CommWorldNext();
		mpi_below = MPI.CommWorldPrev();
		MPI.Printf(stderr, "-- Rank %d | Size=%d | %d -> %d -> %d\n", mpi_rank, mpi_size, mpi_below, mpi_rank, mpi_above );
		
		
		
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
		
		auto_ptr<Workspace> pW0(NULL);
		if( mpi_rank == 0 )
		{
			GhostsInfos no_ghosts( Coord(0,0,0), Coord(0,0,0) );
			GhostsSetup without_ghosts( no_ghosts, no_ghosts );
			pW0.reset( new Workspace( full_layout, without_ghosts, full_region, 1, 1, NULL ) );
		}
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup simulation domain
		//
		////////////////////////////////////////////////////////////////////////
		static const char  *var_names[] = { "u", "v", "Lu", "Lv", "h" };
		static const size_t var_count   = (sizeof(var_names)/sizeof(var_names[0]))/2;
		
		Domain domain( full_layout, full_region, var_count, var_names );
		
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// Initialize fields
		//
		////////////////////////////////////////////////////////////////////////
		{
			Fill::function3 f3( cfunctor3( initV ) );
			Fill::with( f3, domain["v"], domain, domain.X, domain.Y, domain.Z );
		}
		domain["u"].set_all( domain, 0.5 + Real(mpi_rank)/mpi_size );
		
		////////////////////////////////////////////////////////////////////////
		//
		// prepare reactions
		//
		////////////////////////////////////////////////////////////////////////
		Reaction     rxn;
		ODE_Function reaction( &rxn, & Reaction::compute );
		
		
		////////////////////////////////////////////////////////////////////////
		//
		// Run...
		//
		////////////////////////////////////////////////////////////////////////
		const Real dt      = 1e-3;
		Real        t      = 0.0;
		domain["h"].set_all( domain, dt * 0.01 );
		Timings    timings = { 0, 0 };
		
		save_all(0,pW0,domain);
		
		size_t counter = 0;
		for( size_t iter=1; iter <= 6000; ++iter )
		{
			t = (iter-1) * dt;
			domain.cycle(dt, MPI, timings);
			(void) domain.reaction(reaction,t,t+dt);
			if( 0 == (iter%20) )
			{
				MPI.Printf0(stderr,"cycle %5lu: COMM: %7.2f, DIFF: %7.2f        \r", iter, (timings.t_comm*1000)/iter, (timings.t_diff*1000)/iter );
				save_all(++counter,pW0,domain);
			}
		}
		MPI.Printf0(stderr,"\n");
		
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
