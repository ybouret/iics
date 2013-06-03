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
            //std::cerr << "radius=" << radius << std::endl;
            //Shape::Circle(b, center, radius);
            Shape::Rotate(b, numeric<Real>::two_pi * alea<Real>() );
        }
    }
    
    // regularize and broadcast is valid
    validate_bubbles(MPI);
    
    if( !is_valid )
    {
        MPI.Printf(stderr,"Invalid Bubble");
        done = true;
        return;
    }
    broadcast_bubbles(MPI);
    segment();
    P.ldz();
    pressurize_bubbles();
    pressurize_contours();
    compute_gradP(MPI);
    junctions.save_dat( "j" + MPI.CommWorldID + ".dat" );
    junctions.save_inside_of(B, "in" + MPI.CommWorldID + ".dat");
}
