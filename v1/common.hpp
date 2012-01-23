////////////////////////////////////////////////////////////////////////////////
//
// Using yocto amd cliff
// ( Cartesian Lightweight Implementation of Fields Fragmentation )
// with mpi
////////////////////////////////////////////////////////////////////////////////


#ifndef IICS_COMMON_INCLUDED
#define IICS_COMMON_INCLUDED 1

//! 3D settings
#include "yocto/cliff/wksp3d.hpp"
#include "yocto/cliff/laplacian.hpp"
#include "yocto/cliff/fill.hpp"

//! ODE
#include "yocto/math/ode/drvck.hpp"

//! MPI wrappers
#include "yocto/mpi/mpi.hpp"


//! Administrativia
#include "yocto/exception.hpp"

using namespace yocto;
using namespace cliff;
using namespace math;

namespace IICS
{
	////////////////////////////////////////////////////////////////////////////
	//
	// types definition for 3D space
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
	
	////////////////////////////////////////////////////////////////////////////
	//
	// types definition ODE
	//
	////////////////////////////////////////////////////////////////////////////
	
	typedef array<Real>               Variables;    //!< for ODE interface
	typedef ode::field<Real>::type    ODE_Function; //!< Function(dfdt,t,f)
	typedef ode::drvck<Real>::type    ODE_Driver;   //!< Runge-Kutta-Cash-Karp step
	
	
#define IICS_REAL MPI_DOUBLE
	extern int mpi_rank;
	extern int mpi_size;
	extern int mpi_above;
	extern int mpi_below;
	
	struct Timings 
	{
		double t_comm; //!< communication (async+plain)
		double t_diff; //!< diffusion computation
		double t_ode;  //!< ode solving
	};
	
}

#endif
