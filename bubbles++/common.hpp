#ifndef BUBBLES_COMMON_INCLUDED
#define BUBBLES_COMMON_INCLUDED 1

#include "yocto/cliff/wksp3d.hpp"
#include "yocto/cliff/mpi.hpp"
#include "yocto/cliff/fill.hpp"
#include "yocto/wtime.hpp"
#include "yocto/code/rand.hpp"
#include "yocto/auto-ptr.hpp"
#include "yocto/cliff/rwops.hpp"

using namespace yocto;
using namespace cliff;

namespace Bubble
{
	
	typedef double                   Real;          //!< compute in double precision
	typedef fill<Real,Real>          Fill;          //!< to fill arrays
	
	
	extern int mpi_rank;
	extern int mpi_size;
	extern int mpi_last;
	
#define BUBBLE_REAL MPI_DOUBLE
	
}

#endif
