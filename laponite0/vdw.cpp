#include "./vdw.hpp"

namespace Laponite
{
	const Real R_CONSTANT = 8.31;
	
	VanDerWaals:: VanDerWaals( Real _a, Real _b, Real _m, Real _e, Real _t ) throw() :
	A( _a ),
	B( _b ),
	MolarMass( _m ),
	Thickness( _e ),
	Temperature( _t ),
	T_c( 1e2 * 8 * A / ( 27 * R_CONSTANT * B ) ),
	P_c( 1e5 * A/(27*B*B) ),
	V_c( 1e3 * (3*B) ),
	rho_c( MolarMass/V_c ),
	A_Pa( 1e5 * A )
	{
		
	}
	
	VanDerWaals:: ~VanDerWaals() throw() {}
	
	Real VanDerWaals:: pressure( Real surfacic_rho ) const
	{
		const Real rho_SI = surfacic_rho / Thickness;           //!< kg/m^3
		const Real rho    = rho_SI * 1e-3;                      //!< kg/L
		const Real Vm     = MolarMass/rho;                      //!< L/mol
		const Real RT     = R_CONSTANT * Temperature;           //!< J/mol
		if( Vm <= B )
			throw exception("rho=%g => Vm=%g is too high!", surfacic_rho, Vm);
		const Real Vr_SI  = (Vm -B)* 1e-3;                     //!< m^3/mol
		return (RT / Vr_SI) - A_Pa/(Vm*Vm);                    //!< in Pa
	}
	
}