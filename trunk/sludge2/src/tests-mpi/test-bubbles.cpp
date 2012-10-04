#include "yocto/utest/run.hpp"
#include "../bubbles.hpp"

YOCTO_UNIT_TEST_IMPL(bubble)
{
    const mpi & MPI = mpi::init(&argc,&argv);
    const int   rank = MPI.CommWorldRank;
    const int   size = MPI.CommWorldSize;
    
    const PBC pbc(1);
    Bubbles   bubbles( pbc );
    Bubble   *bubble = bubbles.append();
    hashing::sha1 H;
    
    if( MPI.IsMaster)
    {
        const Vertex center( 10*(0.5 - Alea()), 10*(0.5-Alea()));
        const Real   radius = 0.1 + 10 * Alea();
        bubble->map_circle(center, radius);
        bubble->compute_contour();
        bubble->save_dat( "b0.dat" );
    }
    
    MPI.Printf0(stderr,"### Dispatching bubble\n");
    bubble->dispatch(MPI);
    MPI.Printf( stderr, "%d.%d: %08x\n",size,rank, uint32_t(bubble->get_hash(H)));
    MPI.WaitFor(0.1);
    
    
    
    const Real ylo = -5;
    const Real yhi =  5;
    const Real ylen = yhi - ylo;
    const Real ymin = ylo + (rank * ylen)      / size;
    const Real ymax = ylo + ( (rank+1) *ylen ) / size;
    MPI.Printf(stderr, "%d.%d: ymin=%g -> ymax=%g\n", size,rank,ymin,ymax);
    bubble->locate_spots(ymin, ymax);
    MPI.Printf(stderr, "%d.%d: #spots= %u\n", size, rank, unsigned(bubble->spots.size) );
    bubble->save_spots( vformat("s%d.%d.dat",size,rank) );
    for( Spot *spot = bubble->spots.head; spot; spot=spot->next )
    {
        spot->handle->vertex.y += (rank+1) * 0.25;
    }
    
    MPI.Printf0(stderr,"### Assembling bubble\n");
    bubble->assemble( MPI );
    
    if( MPI.IsMaster)
    {
        bubble->save_dat("b1.dat");
    }
}
YOCTO_UNIT_TEST_DONE()

YOCTO_UNIT_TEST_IMPL(bubbles)
{
    const mpi & MPI  = mpi::init(&argc,&argv);
    const int   rank = MPI.CommWorldRank;
    const int   size = MPI.CommWorldSize;
    
    const PBC pbc(10);
    Bubbles   bubbles( pbc );
    hashing::sha1 H;
    
    if( MPI.IsMaster)
    {
        bubbles.create( 2 + 10 * Alea());
        for( Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next )
        {
            const Vertex center( 10*(0.5 - Alea()), 10*(0.5-Alea()));
            const Real   radius = 0.1 + 10 * Alea();
            bubble->map_circle(center, radius);
            bubble->compute_contour();
        }
    }
    
    bubbles.dispatch(MPI);
    
    MPI.Printf( stderr, "%d.%d: %08x\n",size,rank, uint32_t(bubbles.get_hash(H)));
    MPI.WaitFor(0.1);

    const Real ylo  = -5;
    const Real yhi  =  5;
    const Real ylen = yhi - ylo;
    const Real ymin = ylo + (rank * ylen)      / size;
    const Real ymax = ylo + ( (rank+1) *ylen ) / size;
    for( Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next )
    {
        bubble->locate_spots(ymin, ymax);
    }
    
    
    bubbles.assemble(MPI);
    
    
    
}
YOCTO_UNIT_TEST_DONE()