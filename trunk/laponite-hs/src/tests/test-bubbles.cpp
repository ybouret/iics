#include "../bubble.hpp"
#include <iostream>
#include "yocto/string/conv.hpp"
#include "yocto/ios/ocstream.hpp"

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
        assert(b.size>=3);
        
        ios::ocstream fp( vformat("bubble%u.dat", unsigned(b.size)), false );
        const Point *p = b.root;
        for( size_t i=0; i < b.size; ++i, p = p->next )
        {
            fp("%g %g\n", p->x, p->y );
        }
        fp("%g %g\n", p->x, p->y );
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "Unhandled exception" << std::endl;
    }
    return 1;
}