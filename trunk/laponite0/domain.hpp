#ifndef LAPONITE_DOMAIN_INCLUDED 
#define LAPONITE_DOMAIN_INCLUDED 1

#include "./vdw.hpp"

namespace Laponite
{
	
	class Domain : public Workspace
	{
	public:
		static const char  *VarNames[];
		static const size_t VarCount;
		
		const bool   is_first;
		const bool   is_final;
		const Layout bulk; //!< where grad P can be computed
		
		explicit Domain( const Layout &full_layout, const GhostsSetup &setup, const Region &full_region);
		virtual ~Domain() throw();
		
		void exchange_start( const mpi &MPI );
		void exchange_finish( const mpi &MPI );
		void update_boundaries();
		
		//! update rho and compute P
		void compute_P(  const VanDerWaals &gas );
		
	private:
		mpi::Requests requests;
		YOCTO_DISABLE_COPY_AND_ASSIGN(Domain);
	};
}

#endif
