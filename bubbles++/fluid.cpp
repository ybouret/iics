#include "./fluid.hpp"

namespace Bubble
{
	
	Fluid:: ~Fluid() throw() {}
	Fluid::  Fluid() throw() {}
	
	VdW:: VdW( Real temperature ) :
	T( temperature )
	{
		assert( T > 0 );
		
	}

	VdW:: ~VdW() throw() {}
	
	Real VdW:: pressure( Real rho ) const
	{
		static const Real fac = 8.0/3;
		return fac * rho * T / ( 1.0 - rho/3.0 ) - 3*rho*rho;
	}

}