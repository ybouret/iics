#include "../cell.hpp"
#include "yocto/spade/variables.hpp"

void Cell:: dispatch( const mpi &MPI )
{
    MPI.Printf0(stderr, "\tdispatch %u bubbles...\n", unsigned(bubbles.count()));
    bubbles.dispatch(MPI);
    
    MPI.Printf0(stderr, "\tlocate spots...\n");
    for( Bubble *b = bubbles.first(); b;b=b->next)
    {
        b->locate_spots(ymin, ymax);
    }
    
    MPI.Printf0(stderr, "\tsegmentation...\n");
    segmenter.process(bubbles);
    
    MPI.Printf0(stderr, "\tbuild bubble field...\n");
    segmenter.build(B);
    
    MPI.Printf0(stderr, "\tsync bubble field...\n");
    sync1(MPI,B);
    
    variables bvar;
    bvar.append( "B" );
    vtk.save( vformat("b-out-%d.%d.vtk", MPI.CommWorldSize,MPI.CommWorldSize),
             "B full layout", *this,
             bvar,
             outline);
    
    segmenter.save_vtk_n( vformat("j%d%d.vtk",MPI.CommWorldSize,MPI.CommWorldSize), bubbles.lambda);
    
}
