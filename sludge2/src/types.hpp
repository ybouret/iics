#ifndef TYPES_INCLUDED
#define TYPES_INCLUDED 1

#include "yocto/spade/array2d.hpp"
#include "yocto/spade/workspace.hpp"
#include "yocto/spade/rmesh.hpp"
#include "yocto/math/types.hpp"

using namespace yocto;
using namespace spade;
using namespace math;

typedef double                         Real;
typedef complex<Real>                  Complex;
typedef vertex2D<Real>::type           Vertex;
typedef coord2D                        Coord;
typedef layout2D                       Layout;
typedef rmesh<Layout,Real>             Grid;
typedef array1D<Real>                  Array1D;
typedef array2D<Real>                  Array;
typedef array2D<Vertex>                VertexArray;
typedef array1D<Vertex>                VertexArray1D;
typedef ghosts_setup                   GhostsSetup;
typedef fields_setup<Layout>           FieldsSetup;
typedef workspace<Layout,rmesh,Real>   WorkspaceBase;


#if defined(HAS_MPI)
#include "yocto/mpi/mpi.hpp"
#define MPI_REAL_TYPE                  MPI_DOUBLE
#endif


void AleaInit() throw();
Real Alea() throw();

inline Real PBC1( Real x, const Real L, const Real invL ) throw()
{
    static const Real __half = 0.5;
    return x - L * Floor( (invL*x) + __half );
}

//! periodic boundary condition on y
class PBC
{
public:
    const Real L;
    const Real invL;
    const Real lo;
    const Real up;
    PBC( Real length ) throw();
    ~PBC() throw();
    PBC( const PBC &other ) throw();
    
    Real apply( Real y ) const throw();
    
    void operator()( Vertex &v ) const throw(); //!< act on y
    
private:
    YOCTO_DISABLE_ASSIGN(PBC);
};

void SaveGrid( const Grid &, const string &filename );


#endif
