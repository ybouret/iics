#include "yocto/utest/run.hpp"
#include "../mpi/cell.hpp"
#include "yocto/code/utils.hpp"

YOCTO_UNIT_TEST_IMPL(cell)
{
    const mpi &MPI = mpi::init(&argc,&argv);
    const Coord  N(10,20);
    const Vertex L(2.0,3.0);
    
    FieldsSetup F(2);
    Y_SPADE_FIELD(F, "B", Array);
    
    Cell cell(MPI,N,L,F);
    Bubbles &bubbles = cell.bubbles;
    
    //--------------------------------------------------------------------------
    // create the bubble(s) on master
    //--------------------------------------------------------------------------
    if( MPI.IsMaster)
    {
        Bubble *b = bubbles.append();
        Vertex center(L.x/2,0);
        Real   radius = min_of(L.x,L.y)/5;
        b->map_circle(center, radius);
        b->compute_contour();
    }
    
    //--------------------------------------------------------------------------
    //-- broadcast bubbles with data & find spots
    //--------------------------------------------------------------------------
    cell.dispatch(MPI);
    
        
}
YOCTO_UNIT_TEST_DONE()

