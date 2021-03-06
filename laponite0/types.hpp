#ifndef LAPONITE_TYPES_INCLUDED
#define LAPONITE_TYPES_INCLUDED 1

#include "yocto/cliff/wksp2d.hpp"
#include "yocto/cliff/fill.hpp"
#include "yocto/mpi/mpi.hpp"

using namespace yocto;
using namespace cliff;
using namespace math;


namespace Laponite
{
	
	////////////////////////////////////////////////////////////////////////////
	//
	// types definition for 2D space
	//
	////////////////////////////////////////////////////////////////////////////
	typedef double                   Real;         //!< compute in double precision
	typedef coord2D                  Coord;        //!< logical coordinate for indexing
	typedef layout2D                 Layout;       //!< layout for 2D arrays
	typedef region2D<Real>::type     Region;       //!< 2D region
	typedef array2D<Real>            Array;        //!< array of Real, in 2D
	typedef array1D<Real>            Row;          //!< array of Real, in 1D
	typedef vertex2D<Real>::type     Vertex;       //!< vertices or Reals
	typedef wksp2D<Real,Real>        Workspace;    //!< 2D workspace
	typedef ghost<Real,Coord>        Ghost;        //!< ghost of Reals
	typedef fill<Real,Real>          Fill;         //!< to fill arrays
	typedef Fill::function2          FillFunctor;  //!< with this functor
	typedef ghosts_infos<Coord>      GhostsInfos;  //!< count/async
	typedef ghosts_setup<Coord>      GhostsSetup;  //!< infos for lower an upper
	
	
	////////////////////////////////////////////////////////////////////////////
	//
	// for MPI
	//
	////////////////////////////////////////////////////////////////////////////
	#define IICS_REAL MPI_DOUBLE
	extern int mpi_rank;
	extern int mpi_size;
	extern int mpi_prev; 
	extern int mpi_next;
	extern int mpi_last;
}

#endif
