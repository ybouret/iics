#include "yocto/utest/run.hpp"
#include "../bubble.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/code/utils.hpp"

static inline
Real Poly( Real x, const Real xa[], const Real ya[], size_t n )
{
    Real ans = 0;
    for( size_t i=0; i < n; ++i )
    {
        const Real yi = ya[i];
        const Real xi = xa[i];
        Real num = 1;
        Real den = 1;
        for( size_t j=0; j < i; ++j )
        {
            const Real xj = xa[j];
            num *= x  - xj;
            den *= xi - xj;
        }
        for( size_t j=i+1;j<n;++j)
        {
            const Real xj = xa[j];
            num *= x  - xj;
            den *= xi - xj;
        }
        ans += (num*yi)/den;
    }
    return ans;
}

static inline void save_poly( const string &name, const Real xa[], const Real ya[], size_t n )
{
    ios::ocstream fp( name, false );
    const Real   x1 = xa[0];
    const Real   xn = xa[n-1];
    const size_t m  = 100;
    for( size_t i=0; i <=m; ++i)
    {
        const Real x = x1 + i * (xn-x1)  / m;
        fp("%g %g\n",x,Poly(x,xa,ya,n));
    }
    
}

YOCTO_UNIT_TEST_IMPL(poly)
{
    Real xa[4] = { -1, 0, 1, 2 };
    Real ya[4] = {  1, 1, 1, 1 };
    const size_t n = sizeof(xa)/sizeof(xa[0]);
    
    save_poly("poly0.dat",xa,ya,n);
    
    ya[0] = 0;
    save_poly("poly1.dat",xa,ya,n);

    ya[3] = 0;
    save_poly("poly2.dat",xa,ya,n);
    
    ya[3] = 2;
    save_poly("poly3.dat",xa,ya,n);

    
}
YOCTO_UNIT_TEST_DONE()

