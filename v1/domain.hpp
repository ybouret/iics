#ifndef IICS_DOMAIN_INCLUDED
#define IICS_DOMAIN_INCLUDED 1

#include "common.hpp"
#include "yocto/wtime.hpp"

namespace  IICS {
	
	
	class Domain : public Workspace
	{
	public:
		static const size_t VarStart = 1;
		
		static const Coord        GCount;
		static const Coord        GAsync;
		static const GhostsInfos  GInfos;
		static const GhostsSetup  GSetup;
		
		//! constructor
		/**
			\param full_layout layout of the global simulation
			\param full_region region of the global simulation
			\param fields      number of fields
			\param names       2*fields valid names
		 
			'names' must be of the form { "u", "v", ..., "Lu", "Lv", .. }
			The MPI global variables rank, size, above and below must be set.
		 */
		explicit Domain( const Layout &full_layout, 
						 const Region &full_region,
						 const size_t  fields,
						 const char   *names[]
						);
		virtual ~Domain() throw();
		
	
		vector<size_t> field_index; //!< indices of active variables
		vector<size_t> delta_index; //!< indices of corresponding laplacian
		const size_t   num_fields;  //!< field_index.size(), same as delta_index.size()
		mpi::Requests  requests;    //!< 4 * num_fields
		wtime          chrono;
		
		double start_exchanges( const mpi &MPI );
		double start_laplacian( Real dt );
		double finish_exchanges( const mpi &MPI );
		double finish_laplacian( Real dt );
		double update();
		
		void   cycle( Real dt, const mpi &MPI, Timings &timings );
		
		
	private:
		YOCTO_DISABLE_COPY_AND_ASSIGN(Domain);
		
	};
	
}

#endif
