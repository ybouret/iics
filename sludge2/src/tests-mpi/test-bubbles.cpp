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
    }
    
    bubble->dispatch(MPI);
    MPI.Printf( stderr, "%d.%d: %08x\n",size,rank, uint32_t(bubble->get_hash(H)));
    MPI.__WaitFor(0.1);
}
YOCTO_UNIT_TEST_DONE()

YOCTO_UNIT_TEST_IMPL(bubbles)
{
    const mpi & MPI  = mpi::init(&argc,&argv);
    const int   rank = MPI.CommWorldRank;
    const int   size = MPI.CommWorldSize;
    
    const PBC pbc(1);
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
        }
    }
    
    bubbles.dispatch(MPI);
    
    MPI.Printf( stderr, "%d.%d: %08x\n",size,rank, uint32_t(bubbles.get_hash(H)));
    MPI.__WaitFor(0.1);

}
YOCTO_UNIT_TEST_DONE()