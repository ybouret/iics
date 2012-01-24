#ifndef LAPONITE_VDW_INCLUDED
#define LAPONITE_VDW_INCLUDED

#include "./types.hpp"

namespace Laponite
{

	extern const Real R_CONSTANT;
	
	class VanDerWaals
	{
	public:
		explicit VanDerWaals( Real _a, Real _b, Real _m, Real _e, Real _t ) throw();
		virtual ~VanDerWaals() throw();
		
		const Real A;           //! in bar * (L/mol)^2
		const Real B;           //! in L/mol
		const Real MolarMass;   //!< in kg/mol
		const Real Thickness;   //!< in meters
		Real       Temperature; //!< in Kelvin
		
		//! get the pressure in bar
		/**
		 \param surfacic rho in kg/m^2
		 */
		Real pressure( Real surfacic_rho ) const;
		
	private:
		YOCTO_DISABLE_COPY_AND_ASSIGN(VanDerWaals);
	};
	
}


#endif
