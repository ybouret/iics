#include "simulation.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/code/rand.hpp"
#include "../shape.hpp"


void Simulation:: init_one_bubble()
{
    if( MPI.IsFirst )
    {
        bubbles.auto_delete();
        const Real radius = min_of( full_region.length.x, full_region.length.y)/3;
        {
            const Vertex center = full_region.vmin + 0.5 * full_region.length;
            Bubble *b = bubbles.append();
            Shape::Blob(b, center, radius, 0.4 + 0.55 * alea<Real>(), 0.1 + 0.8 * alea<Real>() );
            Shape::Rotate(b, numeric<Real>::two_pi * alea<Real>() );
        }
        bubbles.regularize_all();
        bubbles.head->save_dat("b0.dat");
    }
    broadcast_bubbles(MPI);
    segment();
    reset_pressure();
    junctions.save_dat( "j" + MPI.CommWorldID + ".dat" );
    junctions.save_inside_of(B, "in" + MPI.CommWorldID + ".dat");
}
