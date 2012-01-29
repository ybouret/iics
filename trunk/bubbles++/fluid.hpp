#ifndef BUBBLE_FLUID_INCLUDED
#define BUBBLE_FLUID_INCLUDED 1

#include "./common.hpp"

namespace Bubble
{
	
	class Fluid 
	{
	public:
		virtual ~Fluid() throw();
		virtual Real pressure( Real rho ) const = 0;
		
	protected:
		explicit Fluid() throw();
		
	private:
		YOCTO_DISABLE_COPY_AND_ASSIGN(Fluid);
	};
	
	class VdW : public Fluid
	{
	public:
		Real T;
		virtual Real pressure( Real rho ) const;
		
		explicit VdW(Real temperature);
		virtual ~VdW() throw();
		
	private:
		YOCTO_DISABLE_COPY_AND_ASSIGN(VdW);
	};
	
}

#endif
