#ifndef LAPONITE_DOMAIN_INCLUDED 
#define LAPONITE_DOMAIN_INCLUDED 1

#include "./types.hpp"

namespace Laponite
{
	
	class Domain : public Workspace
	{
	public:
		static const char  *VarNames[];
		static const size_t VarCount;
		
		explicit Domain( const Layout &full_layout, const GhostsSetup &, const Region &full_region);
		virtual ~Domain() throw();
		
		void exchange_start( const mpi &MPI );
		
		
	private:
		mpi::Requests requests;
		YOCTO_DISABLE_COPY_AND_ASSIGN(Domain);
	};
}

#endif
