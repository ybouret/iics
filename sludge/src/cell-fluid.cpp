#include "cell.hpp"

Real   Cell:: Lambda( Real g ) const
{
    return 1;
}

Vertex Cell:: velocity_from( const Vertex &g ) const
{
    return - Lambda( g.norm() ) * g;
}