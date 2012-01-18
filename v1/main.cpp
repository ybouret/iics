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
		
		std::cerr << "-- size = " << MPI.CommWorldSize << std::endl;
		std::cerr << "-- rank = " << MPI.CommWorldRank << std::endl;
		
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
