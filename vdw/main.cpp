#include "yocto/math/fcn/zfind.hpp"
#include "yocto/math/fcn/integrate.hpp"

#include "yocto/exceptions.hpp"
#include "yocto/math/types.hpp"
#include "yocto/code/utils.hpp"

using namespace yocto;
using namespace math;

static inline double pressure_rho( double rho, double T )
{
	return 8.0*rho*T/(3.0 - rho ) - 3.0 *rho*rho;
}

static inline double pressure_V( double V, double T )
{
	return 8.0*T/(3.0 * V - 1.0 ) - 3.0 /(V*V);
}


class Fluid
{
public:
	inline Fluid( double t ) :
	Tr(t),
	rho_spin_G( 2 * ( 1+cos(acos(Tr+Tr-1.0)/3.0 + (2.0*numeric<double>::pi/3.0) ) ) ),
	rho_spin_L( 2 * ( 1+cos(acos(Tr+Tr-1.0)/3.0 + (4.0*numeric<double>::pi/3.0) ) ) ),
	V_spin_G( 1.0 / rho_spin_G ),
	V_spin_L( 1.0 / rho_spin_L ),
	P_spin_min( pressure_rho( rho_spin_L, Tr ) ),
	P_spin_max( pressure_rho( rho_spin_G, Tr ) ),
	zrho( this, & Fluid:: zrho__ ),
	zvol( this, & Fluid:: zvol__ ),
	DeltaP( this, & Fluid:: DeltaP__ ),
	DeltaG( this, & Fluid:: DeltaG__ ),
	Pm(0),
	rho_m_G(0),
	rho_m_L(0)
	{
		if( Tr >= 1 ) throw exception("Invalid Tr=%g", Tr );
	}
	
	inline ~Fluid() throw()
	{
	}
	
	
	const double              Tr;
	const double              rho_spin_G;
	const double              rho_spin_L;
	const double              V_spin_G;
	const double              V_spin_L;
	const double              P_spin_min;
	const double              P_spin_max;
	numeric<double>::function zrho;
	numeric<double>::function zvol;
	numeric<double>::function DeltaP;
	numeric<double>::function DeltaG;
	double                    Pm;
	double                    rho_m_G;
	double                    rho_m_L;
	
	double zrho__( double rho )
	{
		return Pm - pressure_rho( rho, Tr );
	}
	
	double zvol__( double V )
	{
		return Pm - pressure_V(V, Tr );
	}
	
	
	double DeltaP__( double V )
	{
		return Pm - pressure_V( V, Tr );
	}
	
	double DeltaG__( double P )
	{
		Pm = P;
		std::cerr << "Pm=" << Pm << std::endl;
		//-- find gaz maxwell Volumes
		zfind<double> solve = { 1e-7 };
		assert(P>=P_spin_min);
		assert(P<=P_spin_max);
		int    n   = 1;
		double Vlo = 0;
		do
		{
			Vlo = 1.0/(3.0-1.0/(++n));
			//std::cerr << "Vlo = " << Vlo << std::endl;
			//std::cerr << "P   = " << pressure_V( Vlo, Tr ) << " / " << Pm << std::endl;
		} while( Vlo >= V_spin_L || pressure_V( Vlo, Tr ) <= Pm );
		
		const double Vmin = solve( zvol, Vlo, V_spin_L );
		std::cerr << "Vmin=" << Vmin << std::endl;
		
		double Vhi = V_spin_G;
		do
		{
			Vhi += V_spin_G;
			std::cerr << "Vhi = " << Vhi << std::endl;
			std::cerr << "P   = " << pressure_V( Vhi, Tr ) << " / " << Pm << std::endl;
		} while( pressure_V( Vhi, Tr ) >= Pm );
		const double Vmax = solve( zvol, V_spin_G, Vhi );
		std::cerr << "Vmax=" << Vmax << std::endl;
		
		
		//std::cerr << "Vmin= " << Vmin << " -> Vmax= " << Vmax << std::endl;
		
		//integrate<double> intg = { 1e-7 };
		//const double dG = intg( Vmin, Vmax, DeltaP );
		//std::cerr << "dG = " << dG << std::endl;
		//return dG;
		return 0;
	}
	
	void find()
	{
		zfind<double> solve = { 1e-7 };
#if 0
		const double P_range = P_spin_max - P_spin_min;
		const double dG_top  = DeltaG__( P_spin_max );
		
		std::cerr << "dG_top=" << dG_top << std::endl;
		return;
		Pm      =  solve( DeltaG, P_spin_min, P_spin_max );
#endif
	}
	
	
private:
	YOCTO_DISABLE_COPY_AND_ASSIGN(Fluid);
};


int main( int argc, char *argv[] )
{
	
	try
	{
		Fluid fluid(0.7);
		std::cerr << "T          = " << fluid.Tr    << std::endl;
		std::cerr << "rho_spin_L = " << fluid.rho_spin_L << std::endl;
		std::cerr << "rho_spin_G = " << fluid.rho_spin_G << std::endl;

		std::cerr << "V_spin_L   = " << fluid.V_spin_L << std::endl;
		std::cerr << "V_spin_G   = " << fluid.V_spin_G << std::endl;
		
		std::cerr << "P_spin_min = " << fluid.P_spin_min << std::endl;
		std::cerr << "P_spin_max = " << fluid.P_spin_max << std::endl;
		
		fluid.DeltaG( 0.5 * (fluid.P_spin_max + fluid.P_spin_min ) );
							 
		
		fluid.find();
		std::cerr << "P_maxwell  = " << fluid.Pm << std::endl;
		std::cerr << "rho_m_G    = " << fluid.rho_m_G << std::endl;
		std::cerr << "rho_m_L    = " << fluid.rho_m_L << std::endl;
		
	}
	catch(...)
	{
		std::cerr << "**** unhandled exception!" << std::endl;
	}
	return -1;
	
}