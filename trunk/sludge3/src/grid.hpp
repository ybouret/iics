#ifndef SLUDGE_GRID_INCLUDED
#define SLUDGE_GRID_INCLUDED 1

#include "yocto/spade/rmesh.hpp"
#include "yocto/spade/region2d.hpp"
#include "yocto/type-ints.hpp"
#include "bubble.hpp"

using namespace spade;

typedef layout2D Layout;
typedef layout1D Layout1D;
typedef coord2D  Coord;

typedef rmesh<Layout,Real>   Grid;
typedef array1D<Real>        Array1D;
typedef region2D<Real>::type Region2D;

#define SLUDGE_INSIDE (0)
#define SLUDGE_TOP    (1)
#define SLUDGE_BOTTOM (2)
#define SLUDGE_LEFT   (4)
#define SLUDGE_RIGHT  (8)

#define SLUDGE_TOP_OR_BOTTOM ( SLUDGE_TOP  | SLUDGE_BOTTOM )
#define SLUDGE_LEFT_OR_RIGHT ( SLUDGE_LEFT | SLUDGE_RIGHT  )

#define SLUDGE_INVALID_COORD (limit_of<unit_t>::maximum)

struct __Grid
{
    
    //! save to gnuplot format
    static void SaveDat( const Grid &grid, const string &fn );
    
    //! Get lambda for bubble resolution
    static Real ComputeLambda( const Grid &grid) throw();
    
    //! locate point on grid
    /**
     \param p any point
     \param lo lower coordinates such 
     that X[lo.x] <= p.x < X[lo.x+1]
     and  Y[lo.y] <= p.y < Y[lo.y+1]
     \return status
     */
    static int Locate( const Grid &grid, const Vertex &p, Coord &lo) throw();
};

#endif
