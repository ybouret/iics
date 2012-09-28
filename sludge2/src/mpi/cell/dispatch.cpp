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
    segmenter.build_bubbles_field(B);
    
//#define SAVE_B
#if defined(SAVE_B)
    variables bvar;
    bvar.append( "B" );
    vtk.save( vformat("b-core-%d.%d.vtk", MPI.CommWorldSize,MPI.CommWorldRank),
             "B full layout", *this,
             bvar,
             outline);
#endif
    
    MPI.Printf0(stderr, "\tsync bubble field...\n");
    sync1(MPI,B);
    
#if defined(SAVE_B)
    vtk.save( vformat("b-sync-%d.%d.vtk", MPI.CommWorldSize,MPI.CommWorldRank),
             "B full layout", *this,
             bvar,
             outline);
    segmenter.save_vtk_n( vformat("j%d-%d.vtk",MPI.CommWorldSize,MPI.CommWorldRank), bubbles.lambda);
#endif
    
    MPI.Printf0(stderr,"\tbuilding effective pressure...\n");
    segmenter.build_effective_pressure(B, Penter, Pleave);
    sync1(MPI,Penter);
    sync1(MPI,Pleave);
}
