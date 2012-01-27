#ifndef BUBBLES_COMMON_INCLUDED
#define BUBBLES_COMMON_INCLUDED 1

#include "yocto/cliff/wksp3d.hpp"
#include "yocto/cliff/mpi.hpp"
#include "yocto/cliff/fill.hpp"
#include "yocto/wtime.hpp"

using namespace yocto;
using namespace cliff;

namespace Bubble
{
	
	typedef double                   Real;          //!< compute in double precision
	typedef coord3D                  Coord;         //!< logical coordinate for indexing
	typedef layout3D                 Layout;        //!< layout for 3d arrays
	typedef region3D<Real>::type     Region;        //!< 3D region
	typedef array3D<Real>            Array;         //!< array of Real, in 3D
	typedef vertex3D<Real>::type     Vertex;        //!< vertices or Reals
	typedef wksp3D<Real,Real>        Workspace;     //!< 3D workspace
	typedef ghost<Real,Coord>        Ghost;         //!< ghost of Reals
	//typedef laplacian<Real,Real>     Laplacian;     //!< to compute laplacians
	typedef fill<Real,Real>          Fill;          //!< to fill arrays
	typedef ghosts_infos<Coord>      GhostsInfos;   //!< count/async
	typedef ghosts_setup<Coord>      GhostsSetup;   //!< infos for lower an upper
	
	extern int mpi_rank;
	extern int mpi_size;
	extern int mpi_last;
	extern int mpi_above;
	extern int mpi_below;
	
}

#endif
