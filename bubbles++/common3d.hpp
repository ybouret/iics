#ifndef COMMON3D_INCLUDED
#define COMMON3D_INCLUDED 1

#include "./common.hpp"

namespace Bubble
{
	typedef coord3D                  Coord;         //!< logical coordinate for indexing
	typedef layout3D                 Layout;        //!< layout for 3d arrays
	typedef region3D<Real>::type     Region;        //!< 3D region
	typedef array3D<Real>            Array;         //!< array of Real, in 3D
	typedef vertex3D<Real>::type     Vertex;        //!< vertices or Reals
	typedef wksp3D<Real,Real>        Workspace;     //!< 3D workspace
	typedef ghost<Real,Coord>        Ghost;         //!< ghost of Reals
	typedef ghosts_infos<Coord>      GhostsInfos;   //!< count/async
	typedef ghosts_setup<Coord>      GhostsSetup;   //!< infos for lower an upper
	
	extern int mpi_above;
	extern int mpi_below;

}

#endif
