#ifndef SLUDGE_MPI_PARAMETERS_INCLUDED
#define SLUDGE_MPI_PARAMETERS_INCLUDED 1

#include "../junctions.hpp"
#include "./bubbles.hpp"
#include "yocto/spade/dataspace.hpp"

typedef fields_setup<Layout> FieldsSetup;
typedef ghosts_setup         GhostsSetup;
typedef array2D<Vertex>      VertexArray;

class Parameters
{
public:
    static const size_t NumGhosts = 2;
    
    explicit Parameters(const mpi    &MPI,
                        const Coord  &N,
                        const Vertex &Q
                        );
    virtual ~Parameters() throw();
    
    FieldsSetup    F;
    GhostsSetup    G;
    const Layout   full_layout;
    const Layout   sim_layout;
    const Vertex   delta;
    const Vertex   inv_delta;
    const Vertex   order1fac;
    const Region2D full_region;
    const Region2D sim_region;
    const unit_t   bulk_imin;
    const unit_t   bulk_imax;
    const unit_t   bulk_jmin;    //!< take care of side
    const unit_t   bulk_jmax;    //!< take care of side
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};

#endif

