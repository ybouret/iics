#include "yocto/utest/run.hpp"
#include "../parameters.hpp"


YOCTO_UNIT_TEST_IMPL(param)
{
    const mpi &MPI = mpi::init(&argc,&argv);
    const Coord  N(10,20);
    const Vertex L(2.0,3.0);
    Parameters   param(N,L,MPI.CommWorldRank,MPI.CommWorldSize);
    array_db     adb;
    Grid         grid(param.sim_layout,adb);
    
    param.setup_grid(grid);
    
    SaveGrid(grid, vformat("g%d.%d.dat", MPI.CommWorldSize, MPI.CommWorldRank) );
    
    MPI.WaitFor(MPI.CommWorldRank*0.5);
    std::cerr << MPI.CommWorldSize << "." << MPI.CommWorldRank << std::endl;
    std::cerr << "full_layout: " << param.full_layout << std::endl;
    std::cerr << "sim_layout : " << param.sim_layout << std::endl;
    
    std::cerr.flush();
    MPI.Barrier(MPI_COMM_WORLD);
}
YOCTO_UNIT_TEST_DONE()

