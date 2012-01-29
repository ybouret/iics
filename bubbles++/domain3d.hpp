#ifndef BUBBLE_DOMAIN_INCLUDED
#define BUBBLE_DOMAIN_INCLUDED 1

#include "./fluid.hpp"

namespace Bubble 
{
	
	class Domain : public Workspace
	{
	public:
		explicit Domain(const Layout      &L,
						const GhostsSetup &G,
						const Region      &R );
		virtual ~Domain() throw();
		
		static const char  *VarNames[]; // { "rho", "U", "V", "P", ... }
		static const size_t VarCount;   // sizeof(VarNames)/sizeof(VarNames[0])
		
		vector<size_t> field_index; //!< to be transferred
		
		
		double exchanges_start();
		double exchanges_finish();
		double compute_pressure( const Fluid &F); //! not on ghosts
		double update_fields( Real dt );
		
		wtime chrono;
	private:
		mpi::Requests requests;
		YOCTO_DISABLE_COPY_AND_ASSIGN(Domain);
		
	};
	
}

#endif
