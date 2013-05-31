#include "yocto/utest/run.hpp"
#include "mpi/workspace.hpp"
#include "shape.hpp"
#include "yocto/code/rand.hpp"

YOCTO_UNIT_TEST_IMPL(work)
{
    YOCTO_MPI;
    Workspace W(MPI,Coord(20,30),Vertex(5,5));
    
#if 0
    MPI.WaitFor(0.1*MPI.CommWorldRank);
    std::cerr << "X=";
    for(unit_t i=W.X.lower; i <= W.X.upper; ++i)
        std::cerr << " " << W.X[i];
    std::cerr << std::endl;
    
    std::cerr << "Y=";
    for(unit_t i=W.Y.lower; i <= W.Y.upper; ++i)
        std::cerr << " " << W.Y[i];
    std::cerr << std::endl;
#endif
    
    __Grid::SaveDat(W.mesh, vformat("g%d.%d.dat", MPI.CommWorldSize, MPI.CommWorldRank) );
    
    const Vertex center = 0.5 * (W.full_region.vmin+W.full_region.vmax);
    
    if(MPI.IsFirst)
    {
        {
            Bubble *b = W.bubbles.append();
            Shape::Blob(b, center, W.full_region.length.x/2, 0.5 + 0.4 * alea<Real>(), 0.8 * alea<Real>() );
            Shape::Rotate(b, alea<Real>() * numeric<Real>::two_pi);
        }
        
        W.bubbles.regularize_all();
        for(const Bubble *b = W.bubbles.head;b;b=b->next)
        {
            b->save_all( vformat("b%u", unsigned(b->UID)) );
        }
    }
    
    W.broadcast_bubbles(MPI);
    W.segment();
    W.junctions.save_all( "j" + MPI.CommWorldID );
    W.junctions.save_inside_of( W.B, "in" + MPI.CommWorldID + ".dat" );
    
    MPI.Printf(stderr, "--- Done\n");
    MPI.WaitFor(0.2);
}
YOCTO_UNIT_TEST_DONE()

