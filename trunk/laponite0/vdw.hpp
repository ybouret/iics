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
		const Real T_c;         //!< in Kelvin
		const Real P_c;         //!< in Pa
		const Real V_c;         //!< in m^3/mol
		const Real rho_c;       //!< in kg/m^3
		
		//! get the pressure in Pa
		/**
		 \param surfacic rho in kg/m^2
		 */
		Real pressure( Real surfacic_rho ) const;
		
	private:
		const Real A_Pa; // A * 1e5
		YOCTO_DISABLE_COPY_AND_ASSIGN(VanDerWaals);
	};
	
}


#endif
