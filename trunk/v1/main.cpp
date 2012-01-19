#include "common.hpp"

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
		const size_t rank = MPI.CommWorldSize;
		const size_t size = MPI.CommWorldSize;
		
		std::cerr << "-- size = " << size << std::endl;
		std::cerr << "-- rank = " << rank << std::endl;
		
		////////////////////////////////////////////////////////////////////////
		//
		// setup simulation
		//
		////////////////////////////////////////////////////////////////////////
		
		
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
