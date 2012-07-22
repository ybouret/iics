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
        fp("%g %g\n", p->vertex.x, p->vertex.y );
    }
    fp("%g %g\n", p->vertex.x, p->vertex.y );
}

static inline void save_curvature( const Bubble &b, int level )
{
    assert(b.size>=3);
    
    ios::ocstream fp( vformat("curv%d.dat",level), false );
    const Point *p = b.root;
    for( size_t i=0; i < b.size; ++i, p = p->next )
    {
        fp("%u %g\n", unsigned(i), p->kappa );
    }
    
}

static inline void save_spots( const Bubble &b, int level )
{
    assert(b.size>=3);
    
    ios::ocstream fp( vformat("spot%d.dat",level), false );
    const Spot *spot = b.spots.head;
    for( size_t i=0; i < b.spots.size; ++i, spot = spot->next )
    {
        Point *p = spot->point;
        fp("%g %g\n", p->vertex.x, p->vertex.y );
    }
}



int main( int argc, char * argv[] )
{
    AleaInit();
    try 
    {
        double radius = 1;
        if( argc > 1 )
        {
            radius = strconv::to_real<Real>( argv[1], "radius" );
        }
        Point::Pool pcache;
        Spot::Pool  scache;
        GridMarker::Pool gcache;
        
        const V2D   Length(100,100);
        Real        lambda = 1;
        Bubble      b(0,Length,lambda,pcache,scache,gcache);
        
        b.map_circle( V2D(0,0), radius);
        save_bubble(b,0);
        b.mark_and_find_spots_within(-2, 2);
        save_spots(b, 0);
        std::cerr << "Update 1/2" << std::endl;
        
        b.update_contour();
        b.compute_values();
        save_curvature(b, 0);
        std::cerr << "Area=" << b.area << std::endl;
        
        {
            Point *p = b.root;
            for( size_t i=0; i < (3*b.size)/4; ++i,p=p->next )
            {
                p->vertex *= 2.0 + 3.0 * Alea();
            }
        }
        save_bubble(b,1);
        b.mark_and_find_spots_within(-2, 2);
        save_spots(b, 1);
        std::cerr << "Update 2/2" << std::endl;
        b.update_contour();
        b.compute_values();        
        save_curvature(b, 1);
        std::cerr << "Area=" << b.area << std::endl;
        save_bubble(b,2);
        b.mark_and_find_spots_within(-2, 2);
        save_spots(b, 2);
        return 0;
    }
    catch(...)
    {
        std::cerr << "Unhandled exception" << std::endl;
    }
    return 1;
}