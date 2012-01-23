#include "reaction.hpp"

namespace IICS
{
	Reaction:: ~Reaction() throw()
	{
	}
	
	Reaction:: Reaction() :
	alpha(1),
	beta(1),
	gamma(1),
	delta(1)
	{
	}
	
	void Reaction:: compute( Variables &dfdt, Real T, const Variables &f )
	{
		
		
		const Real u = f[1];
		const Real v = f[2];
		
		Real &dudt = dfdt[1];
		Real &dvdt = dfdt[2];
		
		dudt =  u * ( alpha - beta * v);
		dvdt = -v * ( delta - gamma * u );
		
	}
	
}
