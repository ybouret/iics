#include "../bubbles.hpp"

#include <iostream>
#include "yocto/string/conv.hpp"
#include "yocto/ios/ocstream.hpp"


static inline void info_bubbles( const Bubbles &bubbles, int level )
{
    ios::ocstream fp( vformat("b_info%d.dat", level), false );
    fp("#bubbles= %u\n", unsigned(bubbles.count()));
    unsigned i=1;
    for( const Bubble *b = bubbles.first(); b; b=b->next,++i)
    {
        fp("#\t[[ %d ]]: #points=%u\n", i, unsigned(b->size) ); 
        const Point *p = b->root;
        for( size_t j=b->size;j>0;--j,p=p->next )
        {
            fp("%.5f %.5f\n", p->vertex.x, p->vertex.y);
        }
        fp("%.5f %.5f\n", p->vertex.x, p->vertex.y);
    }
    
}

int main( int argc, char *argv[] )
{
    
    try 
    {
        mpi &MPI = mpi::init( &argc, &argv);
        double radius = 1;
        if( argc > 1 )
        {
            radius = strconv::to_real<Real>( argv[1], "radius" );
        }
        const V2D length(100,100);
        Bubbles bubbles(length);
        if( MPI.IsMaster )
        {
            for( size_t i=5;i>0;--i)
            {
                Bubble *b = bubbles.append();
                b->map_circle( V2D(0.5-Alea(),0.5-Alea()), 1 + 4*Alea() );
            }
        }
        MPI.Printf0( stderr, "Dispatching all %u bubbles\n", unsigned( bubbles.count() ) );
        bubbles.dispatch_all(MPI);
        MPI.Printf(stderr, "Rank %d> #bubbles=%u\n", MPI.CommWorldRank, unsigned( bubbles.count() ) );
        
        info_bubbles(bubbles, MPI.CommWorldRank);
        const double y0 = -6;
        const double H  = 12;
        const double dH = H / MPI.CommWorldSize;
        const double y_lo = y0 + MPI.CommWorldRank * dH;
        const double y_hi = y_lo + dH;

        
        //-- do something with bubbles
        size_t total_count = 0;
        for( Bubble *b = bubbles.first(); b; b=b->next)
        {
            b->mark_and_find_spots_within(y_lo, y_hi);
            total_count += b->spots.size;
        }
        MPI.Printf(stderr, "Rank %d> total #spots=%u\n", MPI.CommWorldRank, unsigned(total_count) );
        
        //-- and collect changes
        bubbles.assemble_all(MPI);
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "Unhandled exception" << std::endl;
    }
    return 1;
}
