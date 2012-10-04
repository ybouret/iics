#include "yocto/utest/run.hpp"
#include "../mpi/cell.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/auto-ptr.hpp"
#include "yocto/spade/mpi/collect0.hpp"
#include "yocto/spade/vtk/writer.hpp"
#include "yocto/spade/variables.hpp"

YOCTO_UNIT_TEST_IMPL(cell)
{
    const mpi &MPI = mpi::init(&argc,&argv);
    const Coord  N(30,40);
    const Vertex L(3.0,4.0);
   
    Cell     cell(MPI,N,L);
    Bubbles &bubbles = cell.bubbles;
    auto_ptr<Workspace> pW;
    Array              *pB = 0;
    Array              *pP = 0;
    vtk_writer          vtk;
    variables           var;
    
    //--------------------------------------------------------------------------
    // create the bubble(s) on master
    //--------------------------------------------------------------------------
    if( MPI.IsFirst)
    {
        //pA.reset( new standalone<Array>( cell.full_layout));
        var.append("B");
        var.append("P");
        FieldsSetup FA(2);
        Y_SPADE_FIELD(FA, "B", Array);
        Y_SPADE_FIELD(FA, "P", Array);
        const GhostsSetup no_ghosts;
        pW.reset( new Workspace(cell.full_layout,FA,no_ghosts) );
        cell.setup_grid( pW->mesh );
        pB = & ((*pW)["B"].as<Array>());
        pP = & ((*pW)["P"].as<Array>());
        Bubble *b      = bubbles.append();
        Vertex  center(L.x/2,0);
        Real    radius = min_of(L.x,L.y)/4;
        b->map_circle(center, radius);
        b->set_pressure(1.0);
        b->compute_contour();
    }
    
    //--------------------------------------------------------------------------
    //-- broadcast bubbles with data & find spots
    //--------------------------------------------------------------------------
    cell.dispatch(MPI);
    cell.compute_pressure(MPI);
    mpi_collect0::get(MPI, pB, cell.B, cell.full_layout);
    mpi_collect0::get(MPI, pP, cell.P, cell.full_layout);
    if( MPI.IsFirst)
    {
        vtk.save("org.vtk", "org", *pW, var, cell.full_layout);
    }
    
    
    SaveGrid( cell.mesh, vformat("g%d.%d.dat", MPI.CommWorldSize,MPI.CommWorldRank));
    
    cell.bubbles.first()->save_dat( vformat("b%d.%d.dat", MPI.CommWorldSize,MPI.CommWorldRank));
    cell.segmenter.save( vformat("j%d.%d.dat", MPI.CommWorldSize,MPI.CommWorldRank));    
    cell.save_B( vformat("i%d.%d.dat", MPI.CommWorldSize,MPI.CommWorldRank));
    cell.segmenter.save_vtk_n( vformat("j%d.%d.vtk", MPI.CommWorldSize,MPI.CommWorldRank), bubbles.lambda/2);
    cell.bubbles.first()->save_vtk( vformat("b%d.%d.vtk", MPI.CommWorldSize,MPI.CommWorldRank));
    cell.bubbles.first()->save_vtk_n( vformat("b%d.%d-n.vtk", MPI.CommWorldSize,MPI.CommWorldRank));
    cell.bubbles.first()->save_vtk_t( vformat("b%d.%d-t.vtk", MPI.CommWorldSize,MPI.CommWorldRank));

    if( MPI.IsFirst)
    {
        Bubble *b      = bubbles.first();
        Vertex  center(L.x/2,L.y/2);
        Real    radius = min_of(L.x,L.y)/5;
        b->map_circle(center, radius);
        b->set_pressure(1.0);
        b->compute_contour();
    }
    
    cell.dispatch(MPI);
    cell.P.ldz();
    cell.compute_pressure(MPI);
    mpi_collect0::get(MPI, pB, cell.B, cell.full_layout);
    mpi_collect0::get(MPI, pP, cell.P, cell.full_layout);
    if( MPI.IsFirst)
    {
        vtk.save("mov.vtk", "mov", *pW, var, cell.full_layout);
    }
    cell.bubbles.first()->save_dat( vformat("n-b%d.%d.dat", MPI.CommWorldSize,MPI.CommWorldRank));
    cell.segmenter.save( vformat("n-j%d.%d.dat", MPI.CommWorldSize,MPI.CommWorldRank));
    cell.save_B( vformat("n-i%d.%d.dat", MPI.CommWorldSize,MPI.CommWorldRank));

        
}
YOCTO_UNIT_TEST_DONE()

