#ifndef IICS_REACTION_INCLUDED
#define IICS_REACTION_INCLUDED 1

#include "common.hpp"

namespace IICS
{

	class Reaction
	{
	public:
		explicit Reaction();
		virtual ~Reaction() throw();
		
		void compute( Variables &dfdt, Real T, const Variables &f );
		
		Real alpha, beta, gamma, delta;
		
	private:
		YOCTO_DISABLE_COPY_AND_ASSIGN(Reaction);
	};
	
	
}

#endif

