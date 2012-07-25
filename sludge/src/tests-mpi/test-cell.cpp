#include "yocto/utest/run.hpp"
#include "../cell.hpp"

#include "yocto/ios/ocstream.hpp"

YOCTO_UNIT_TEST_IMPL(cell)
{
    
    const mpi &MPI = mpi::init( &argc, &argv);
        
    Vertex box(10,10);
    Vertex center(box.x/2,0);
    
    Cell cell(20,30,box,MPI);
    

}
YOCTO_UNIT_TEST_DONE()
