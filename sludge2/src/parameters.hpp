#ifndef PARAMETERS_INCLUDED
#define PARAMETERS_INCLUDED 1

#include "segmenter.hpp"

class Parameters
{
public:
    explicit Parameters(const Coord  &N,
                        const Vertex &L,
                        int           rank,
                        int           size
                        );
    virtual ~Parameters() throw();
    const Layout full_layout; //!< -NY/2 -> NY/2-1
    const Vertex full_length; //!< Lx,Ly
    const Vertex delta;       //!< dX = Lx/full_layout.width;
    const Vertex inv_delta;   //!< 1/delta
    const Vertex inv_delsq;   //!< 1/delta^2
    const Real   rb_factor;  //!< 1.0/(-2*inv_delsq.x - 2*inv_delsq.y): Red/Black factor
    const Vertex inv_twodel;  //!< 1/(2*delta)
    const PBC    pbc;         //!< from LY
    const Layout sim_layout;  //!< from rank and size
    ghosts_setup sim_ghosts;  //!< from rank and size
    FieldsSetup  sim_fields;  //!< to be registered
    void setup_grid( Grid &grid ) const;
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};

#endif

