#include "parameters.hpp"

Parameters:: ~Parameters() throw() {}


static unit_t __half_NY( unit_t NY )
{
    while( NY & 1) ++NY;
    return NY >> 1;
}

static inline
Coord __compute_lo(unit_t NY)
{
    return Coord(0,-__half_NY(NY) );
}

static inline
Coord __compute_up(unit_t NX, unit_t NY)
{
    return Coord(NX,__half_NY(NY)-1);
}

Parameters:: Parameters(const Coord  &N,
                        const Vertex &L,
                        int           rank,
                        int           size) :
full_layout( __compute_lo(N.y), __compute_up(N.x,N.y) ),
full_length( L.x, L.y ),
delta( full_length.x/(full_layout.width.x-1), full_length.y/(full_layout.width.y) ),
inv_delta( 1.0 / delta.x, 1.0/delta.y ),
inv_delsq( inv_delta.x * inv_delta.x, inv_delta.y * inv_delta.y),
rb_factor( 1.0/ (-2*inv_delsq.x -2*inv_delsq.y)),
inv_twodel( 0.5/delta.x, 0.5/delta.y),
pbc(full_length.y),
sim_layout( full_layout.split(rank, size) ),
sim_ghosts(),
sim_fields(4*sizeof(Real))
{
    assert(size>0);
#if 0
    std::cerr << "full_layout: " << full_layout << std::endl;
    std::cerr << "full_length: " << full_length << std::endl;
    std::cerr << "delta      : " << delta << std::endl;
    std::cerr << "lower y    : " << delta.y * full_layout.lower.y << std::endl;
    std::cerr << "upper y    : " << delta.y * full_layout.upper.y << std::endl;
    std::cerr << "upper y+1  : " << delta.y * (full_layout.upper.y+1) << std::endl;
    std::cerr << "PBC        : " << pbc.lo << " -> " << pbc.up << std::endl;
    std::cerr << "sim_layout : " << sim_layout << std::endl;
#endif
    
    //--------------------------------------------------------------------------
    // ghosts layout
    //--------------------------------------------------------------------------
    if( size > 1 )
    {
        //-- parallel
        const int last = size-1;
        sim_ghosts.set_async(ghost::at_lower_y, 2, rank <= 0    ?  last : (rank-1) );
        sim_ghosts.set_async(ghost::at_upper_y, 2, rank >= last ?  0    : (rank+1) );
    }
    else
    {
        //-- not parallel
        sim_ghosts.set_local(on_y, 2);
    }
    
    //--------------------------------------------------------------------------
    // fields
    //--------------------------------------------------------------------------
    FieldsSetup &F = sim_fields;
    Y_SPADE_FIELD(F, "P",     Array);
    Y_SPADE_FIELD(F, "B",     Array);
    Y_SPADE_FIELD(F, "gradP", VertexArray);
    Y_SPADE_FIELD(F, "U",     VertexArray);

}


void Parameters:: setup_grid( Grid &grid ) const
{
    Array1D &X = grid.X();
    //--------------------------------------------------------------------------
    //-- generic X part
    //--------------------------------------------------------------------------
    for( unit_t i=X.lower; i <= X.upper; ++i )
    {
        X[i] = i * delta.x;
    }
    
    //-- specific X part
    if( X.has( full_layout.upper.x) )
        X[ full_layout.upper.x ] = full_length.x;
    
    Array1D &Y = grid.Y();
    //--------------------------------------------------------------------------
    //-- generic Y part
    //--------------------------------------------------------------------------
    for( unit_t j= Y.lower; j <=  Y.upper; ++j )
    {
        Y[j] = j * delta.y;
    }
    
    //-- specific Y part
    if( Y.has(full_layout.lower.y) )
    {
        Y[full_layout.lower.y] = pbc.lo;
    }
    
    if( Y.has(full_layout.upper.y+1) )
    {
        Y[full_layout.upper.y+1] = pbc.up;
    }
    
}


