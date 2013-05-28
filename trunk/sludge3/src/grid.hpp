#ifndef SLUDGE_GRID_INCLUDED
#define SLUDGE_GRID_INCLUDED 1

#include "yocto/spade/rmesh.hpp"
#include "yocto/spade/region2d.hpp"
#include "bubble.hpp"

using namespace spade;

typedef layout2D Layout;
typedef layout1D Layout1D;
typedef coord2D  Coord;

typedef rmesh<Layout,Real>   Grid;
typedef array1D<Real>        Array1D;
typedef region2D<Real>::type Region2D;


struct __Grid
{
    static void SaveDat( const Grid &grid, const string &fn );
};

#endif
