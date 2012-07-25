#ifndef PARAMETERS_INCLUDED
#define PARAMETERS_INCLUDED 1

#include "types.hpp"
#include "yocto/mpi/mpi.hpp"

//! base class for simulation
class Parameters : public FieldsSetup
{
public:
    //! declare all fields
    explicit Parameters(unit_t      Nx, 
                        unit_t      Ny,
                        const Vertex &box,
                        const mpi    &mpi_ref
                        );
    
    //! cleanup
    virtual ~Parameters() throw();
    
    const mpi   &MPI;
    const int    sim_rank;
    const int    sim_size;
    const bool   sim_master;
    const bool   sim_parallel;   
    const Coord  sim_lower;        //!< 0,0
    const Coord  sim_upper;        //!< Nx, Ny
    const Vertex sim_box;          //!< Lx,Ly
    const Layout sim_layout;       //!< Lower,Upper
    const Vertex sim_lower_corner; //!< (0,-Ly/2)
    const Vertex sim_upper_corner; //!< (Lx,Ly/2)
    const Region sim_region;       //!< BotLeft->TopRight
    const Layout sub_layout;       //!< splitted
    const Region sub_region;       //!< splited from SubLayout
    GhostsSetup  gs;               //!< info about ghosts
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};




#endif
