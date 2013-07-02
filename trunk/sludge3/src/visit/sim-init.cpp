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
        
        Real radius = min_of( full_region.length.x, full_region.length.y)/3;
        {
            const Vertex center = full_region.vmin + 0.5 * full_region.length;
            
            
            if(flag == "cc")
            {
                Shape::Circle(bubbles.append(), center, radius/2);
                goto PREPARE;
            }
            
            if(flag == "sq")
            {
                Shape::Square(bubbles.append(), center, radius);
                goto PREPARE;
            }
            
            if( flag == "m" )
            {
                radius /= 2;
                {
                    const Vertex C = full_region.vmin + 0.25 * full_region.length;
                    Shape::Blob(bubbles.append(), C, radius, 0.4 + 0.55 * alea<Real>(), 0.1 + 0.8 * alea<Real>() );
                }
                
                {
                    const Vertex C = full_region.vmin + 0.75 * full_region.length;
                    Shape::Blob(bubbles.append(), C, radius, 0.4 + 0.55 * alea<Real>(), 0.1 + 0.8 * alea<Real>() );
                }
                
                {
                    const Vertex C = full_region.vmin + Vertex(0.25*full_region.length.x,0.75*full_region.length.y);
                    Shape::Blob(bubbles.append(), C, radius, 0.4 + 0.55 * alea<Real>(), 0.1 + 0.8 * alea<Real>() );
                }
                
                {
                    const Vertex C = full_region.vmin + Vertex(0.75*full_region.length.x,0.25*full_region.length.y);
                    Shape::Blob(bubbles.append(), C, radius, 0.4 + 0.55 * alea<Real>(), 0.1 + 0.8 * alea<Real>() );
                }
                
                goto PREPARE;
            }
            
            // default
            Shape::Blob(bubbles.append(), center, radius, 0.4 + 0.55 * alea<Real>(), 0.1 + 0.8 * alea<Real>() );
            Shape::Rotate(bubbles.head, numeric<Real>::two_pi * alea<Real>() );
        }
    }
    
PREPARE:;
    
#if 0
    junctions.save_dat( "j" + MPI.CommWorldID + ".dat" );
    junctions.save_inside_of(B, "in" + MPI.CommWorldID + ".dat");
    if(MPI.IsFirst)
        bubbles.head->save_all("b0");
#endif
    
}
