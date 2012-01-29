#ifndef COMMON2D_INCLUDED
#define COMMON2D_INCLUDED 1

#include "./common.hpp"

namespace Bubble
{
	typedef coord2D                  Coord;         //!< logical coordinate for indexing
	typedef layout2D                 Layout;        //!< layout for 2D arrays
	typedef region2D<Real>::type     Region;        //!< 2D region
	typedef array2D<Real>            Array;         //!< array of Real, in 2D
	typedef vertex2D<Real>::type     Vertex;        //!< vertices or Reals
	typedef wksp2D<Real,Real>        Workspace;     //!< 2D workspace
	typedef ghost<Real,Coord>        Ghost;         //!< ghost of Reals
	typedef ghosts_infos<Coord>      GhostsInfos;   //!< count/async
	typedef ghosts_setup<Coord>      GhostsSetup;   //!< infos for lower an upper
	
	extern int mpi_prev;
	extern int mpi_next;

}

#endif
