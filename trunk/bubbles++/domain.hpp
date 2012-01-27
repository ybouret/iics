#ifndef BUBBLE_DOMAIN_INCLUDED
#define BUBBLE_DOMAIN_INCLUDED 1

#include "./common.hpp"

namespace Bubble 
{
	
	class Domain : public Workspace
	{
	public:
		explicit Domain(const Layout      &L,
						const GhostsSetup &G,
						const Region      &R );
		virtual ~Domain() throw();
		
		static const char  *VarNames[]; // { "rho", "U", "V",... }
		static const size_t VarCount;   // sizeof(VarNames)/sizeof(VarNames[0])
		
		vector<size_t> field_index; //!< to be transferred
		
		
		double exchanges_start();
		double exchanges_finish();
		
		wtime chrono;
	private:
		mpi::Requests requests;
		YOCTO_DISABLE_COPY_AND_ASSIGN(Domain);
		
	};
	
}

#endif
