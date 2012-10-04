#include "cell.hpp"
#include "yocto/code/utils.hpp"

Cell:: ~Cell() throw()
{
    
}

Cell:: Cell(const mpi         &MPI,
            const Coord       &N,
            const Vertex      &L) :
Parameters(N,L,MPI.CommWorldRank,MPI.CommWorldSize),
Workspace( sim_layout, sim_fields, sim_ghosts ),
X( mesh.X() ),
Y( mesh.Y() ),
B(       (*this)["B"     ].as<Array>()       ),
P(       (*this)["P"     ].as<Array>()       ),
gradP(   (*this)["gradP" ].as<VertexArray>() ),
U(       (*this)["U"     ].as<VertexArray>() ),
Penter(  (*this)["Penter"].as<VertexArray>() ),
Pleave(  (*this)["Pleave"].as<VertexArray>() ),
segmenter( mesh ),
bubbles( pbc ),
ymin(0),
ymax(0),
vtk()
{
    
    //--------------------------------------------------------------------------
    //-- compute the grid
    //--------------------------------------------------------------------------
    setup_grid( mesh );
    
    //--------------------------------------------------------------------------
    //-- tune the data
    //--------------------------------------------------------------------------
    (Real&)ymin = Y[ this->lower.y-1 ];
    (Real&)ymax = Y[ this->upper.y+1 ];
    MPI.PrintfI(stderr, "ymin=%g / ymax=%g / pbc=[%g;%g]\n", ymin, ymax, pbc.lo, pbc.up);
    
    //--------------------------------------------------------------------------
    //-- create segments
    //--------------------------------------------------------------------------
    segmenter.create();
    
    //--------------------------------------------------------------------------
    //-- tune bubbles
    //--------------------------------------------------------------------------
    bubbles.lambda = min_of<Real>( delta.x, delta.y )/2;
    
}


void Cell:: save_B( const string &filename ) const
{
    ios::ocstream fp( filename, false);
    for( unit_t j=sim_layout.lower.y;j<=sim_layout.upper.y;++j)
    {
        for(unit_t i=sim_layout.lower.x;i<=sim_layout.upper.x;++i)
        {
            if( B[j][i]>0 )
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
    }
}

void Cell:: save_outB( const string &filename ) const
{
    ios::ocstream fp( filename, false);
    for( unit_t j=Y.lower;j<=Y.upper;++j)
    {
        for(unit_t i=X.lower;i<=X.upper;++i)
        {
            if( B[j][i]>0 )
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
    }
}


void Cell:: legalize(const yocto::mpi &MPI)
{
    if( MPI.IsFirst )
    {
        bubbles.update_topology();
    }
    dispatch(MPI);
    compute_pressure(MPI);
}
