#include "../bubble.hpp"
#include <iostream>
#include "yocto/string/conv.hpp"
#include "yocto/ios/ocstream.hpp"


static inline void save_bubble( const Bubble &b )
{
    assert(b.size>=3);
    
    ios::ocstream fp( vformat("bubble%u.dat", unsigned(b.size)), false );
    const Point *p = b.root;
    for( size_t i=0; i < b.size; ++i, p = p->next )
    {
        fp("%g %g\n", p->x, p->y );
    }
    fp("%g %g\n", p->x, p->y );
}


int main( int argc, char * argv[] )
{
    try 
    {
        double radius = 1;
        if( argc > 1 )
        {
            radius = strconv::to_real<Real>( argv[1], "radius" );
        }
        PCache cache;
        Bubble b(cache);
        
        b.map_circle( V2D(0,0), radius);
        save_bubble(b);
        std::cerr << "Update 1/2" << std::endl;
        
        b.update();
        
        {
            Point *p = b.root;
            for( size_t i=0; i < b.size/2; ++i,p=p->next )
            {
                V2D &v = *p;
                v *= 2.0;
            }
        }
        std::cerr << "Update 2/2" << std::endl;
        b.update();
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "Unhandled exception" << std::endl;
    }
    return 1;
}