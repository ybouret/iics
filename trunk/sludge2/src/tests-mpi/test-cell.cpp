#include "yocto/utest/run.hpp"
#include "../mpi/cell.hpp"

YOCTO_UNIT_TEST_IMPL(cell)
{
    const mpi &MPI = mpi::init(&argc,&argv);
    const Coord  N(10,20);
    const Vertex L(2.0,3.0);
    
    FieldsSetup F(2);
    Y_SPADE_FIELD(F, "B", Array);
    
    Cell cell(MPI,N,L,F);
    
}
YOCTO_UNIT_TEST_DONE()

