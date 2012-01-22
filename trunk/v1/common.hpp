////////////////////////////////////////////////////////////////////////////////
//
// Using yocto amd cliff
// ( Cartesian Lightweight Implementation of Fields Fragmentation )
// with mpi
////////////////////////////////////////////////////////////////////////////////


#ifndef IICS_COMMON_INCLUDED
#define IICS_COMMON_INCLUDED 1

#include "yocto/cliff/wksp3d.hpp"
#include "yocto/cliff/laplacian.hpp"
#include "yocto/cliff/fill.hpp"

#include "yocto/exception.hpp"
#include "yocto/mpi/mpi.hpp"

using namespace yocto;
using namespace cliff;


namespace IICS
{
	////////////////////////////////////////////////////////////////////////////
	//
	// types definition
	//
	////////////////////////////////////////////////////////////////////////////
	typedef double                   Real;         //!< compute in double precision
	typedef coord3D                  Coord;        //!< logical coordinate for indexing
	typedef layout3D                 Layout;       //!< layout for 3d arrays
	typedef region3D<Real>::type     Region;       //!< 3D region
	typedef array3D<Real>            Array;        //!< array of Real, in 3D
	typedef vertex3D<Real>::type     Vertex;       //!< vertices or Reals
	typedef wksp3D<Real,Real>        Workspace;    //!< 3D workspace
	typedef ghost<Real,Coord>        Ghost;        //!< ghost of Reals
	typedef laplacian<Real,Real>     Laplacian;    //!< to compute laplacians
	typedef fill<Real,Real>          Fill;         //!< to fill arrays
	typedef Fill::function3          FillFunctor;  //!< with this functor
	typedef ghosts_infos<Coord>      GhostsInfos;  //!< count/async
	typedef ghosts_setup<Coord>      GhostsSetup;  //!< infos for lower an upper
	
#define IICS_REAL MPI_DOUBLE
	extern int mpi_rank;
	extern int mpi_size;
	extern int mpi_above;
	extern int mpi_below;
	
	struct Timings 
	{
		double t_comm; //!< communication (async+plain)
		double t_diff; //!< diffusion computation
	};
	
}

#endif
