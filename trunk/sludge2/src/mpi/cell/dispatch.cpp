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
    

    
    MPI.Printf0(stderr, "\tsync bubble field...\n");
    sync1(MPI,B);

    
    MPI.Printf0(stderr,"\tbuilding effective pressure...\n");
    segmenter.build_effective_pressure(B, P, Penter, Pleave);
    
    sync1(MPI,Penter);
    sync1(MPI,Pleave);
    
}
