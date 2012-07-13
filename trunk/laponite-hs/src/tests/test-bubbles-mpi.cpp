#include "../bubble.hpp"
#include <iostream>
#include "yocto/string/conv.hpp"
#include "yocto/ios/ocstream.hpp"

static inline void save_bubble( const Bubble &b, int level )
{
    assert(b.size>=3);
    
    ios::ocstream fp( vformat("bubble%d.dat",level), false );
    const Point *p = b.root;
    for( size_t i=0; i < b.size; ++i, p = p->next )
    {
        fp("%g %g\n", p->x, p->y );
    }
    fp("%g %g\n", p->x, p->y );
}

static inline void save_spots( const Bubble &b, int level )
{
    assert(b.size>=3);
    
    ios::ocstream fp( vformat("spot%d.dat",level), false );
    const Spot *spot = b.spots.head;
    for( size_t i=0; i < b.spots.size; ++i, spot = spot->next )
    {
        Point *p = spot->point;
        fp("%g %g\n", p->x, p->y );
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
        Point::Pool pcache;
        Spot::Pool  scache;
        Bubble b(pcache,scache);
        if( MPI.IsMaster )
        {
            b.map_circle( V2D(0,0), radius);
            b.update_points();
        }
        b.dispatch(MPI);
        save_bubble( b, MPI.CommWorldRank );
        
        const double R    = radius * 1.1;
        const double H    = R+R;
        const double dH   = H / MPI.CommWorldSize;
        const double y0   = -R;
        const double y_lo = y0 + MPI.CommWorldRank * dH;
        const double y_hi = y_lo + dH;
        
        b.build_spots(y_lo,y_hi);
        MPI.Printf( stderr, "Rank %d> [%g,%g](+%g): #spots= %u\n", MPI.CommWorldRank, y_lo,y_hi, dH, unsigned(b.spots.size) );
        save_spots(b,MPI.CommWorldRank);
        if( b.active )
        {
            b.update_area();
            //-- move concerned points
            for( Spot *sp = b.spots.head; sp; sp = sp->next )
            {
                Point *p = sp->point;
                p->x += 0.2 * radius * (0.5 - Alea() );
                p->y += 0.2 * radius * (0.5 - Alea() );
            }
        }
        
        b.collect(MPI);
        
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "Unhandled exception" << std::endl;
    }
    return 1;
}