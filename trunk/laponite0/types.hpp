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
	typedef layout2D                 Layout;       //!< layout for 3d arrays
	typedef region2D<Real>::type     Region;       //!< 3D region
	typedef array2D<Real>            Array;        //!< array of Real, in 3D
	typedef vertex2D<Real>::type     Vertex;       //!< vertices or Reals
	typedef wksp2D<Real,Real>        Workspace;    //!< 3D workspace
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
	
	extern int mpi_rank;
	extern int mpi_size;
	extern int mpi_prev; 
	extern int mpi_next;
	
}

#endif
