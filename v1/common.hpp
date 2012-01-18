////////////////////////////////////////////////////////////////////////////////
//
// Using yocto amd cliff
// ( Cartesian Lightweight Implementation of Fields Fragmentation )
// with mpi
////////////////////////////////////////////////////////////////////////////////


#ifndef IICS_COMMON_INCLUDED
#define IICS_COMMON_INCLUDED 1

#include "yocto/cliff/wksp3d.hpp"
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
	typedef double               Real;      //!< compute in double precision
	typedef coord3D              Coord;     //!< logical coordinate for indexing
	typedef layout3D             Layout;    //!< layout for 3d arrays
	typedef array3D<Real>        Array;     //!< array of Real, in 3D
	typedef vertex3D<Real>::type Vertex;    //!< vertices or Reals
	typedef wksp3D<Real,Real>    Workspace; //!< 3D workspace
}

#endif
