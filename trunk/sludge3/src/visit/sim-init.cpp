#include "simulation.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/code/rand.hpp"
#include "../shape.hpp"


void Simulation:: init_one_bubble( const char *kind )
{
    if( MPI.IsFirst )
    {
        bubbles.auto_delete();
        const string flag(kind);
        
        const Real radius = min_of( full_region.length.x, full_region.length.y)/3;
        {
            const Vertex center = full_region.vmin + 0.5 * full_region.length;
            Bubble *b = bubbles.append();
            
            
            if(flag == "cc")
            {
                Shape::Circle(b, center, radius);
                goto PREPARE;
            }
            
            if(flag == "sq")
            {
                 Shape::Square(b, center, radius);
                goto PREPARE;
            }
            
            // default
            Shape::Blob(b, center, radius, 0.4 + 0.55 * alea<Real>(), 0.1 + 0.8 * alea<Real>() );
            Shape::Rotate(b, numeric<Real>::two_pi * alea<Real>() );
        }
    }
    
PREPARE:
    // regularize and broadcast is_valid
    validate_bubbles(MPI);
    
    if( !is_valid )
    {
        MPI.Printf(stderr,"Invalid Bubble");
        done = true;
        return;
    }
    broadcast_bubbles(MPI);
    segment();
    
    junctions.save_dat( "j" + MPI.CommWorldID + ".dat" );
    junctions.save_inside_of(B, "in" + MPI.CommWorldID + ".dat");
    if(MPI.IsFirst)
        bubbles.head->save_all("b0");
    P.ldz();
    pressurize_bubbles();
    pressurize_contours();
    compute_gradP(MPI);
}
