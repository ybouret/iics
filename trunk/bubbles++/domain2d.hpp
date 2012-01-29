#ifndef BUBBLE_DOMAIN2D_INCLUDED
#define BUBBLE_DOMAIN2D_INCLUDED 1

#include "./fluid.hpp"
#include "./common2d.hpp"

namespace Bubble 
{
	
	class Domain : public Workspace
	{
	public:
		static const char * VarNames[];
		static const size_t VarCount;
		static const size_t VarMPI;
		
		explicit Domain(const Layout      &L,
						const GhostsSetup &G,
						const Region      &R );
		virtual ~Domain() throw();
		
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
