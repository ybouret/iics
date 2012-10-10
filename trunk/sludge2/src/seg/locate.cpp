#include "../segmenter.hpp"

void Segmenter:: locate_value( const Real z, const Array1D &Z, unit_t &klo, unit_t &kup ) throw()
{
#if 0
    if(!(z>=Z[Z.lower]))
    {
        fprintf( stderr, "error: ! z=%g >= Z[%ld]=%g\n", z, Z.lower, Z[Z.lower] );
    }
#endif
    
    assert(z>=Z[Z.lower]);
    assert(z<Z[Z.upper]);
    klo = Z.lower;
    kup = Z.upper;
    while( kup - klo > 1 )
    {
        const unit_t k = (kup+klo)>>1;
        if (z<Z[k])
            kup=k;
        else
            klo=k;
    }
    assert(z<Z[kup]);
    assert(z>=Z[klo]);
}

void Segmenter:: locate_vertex( const Vertex &v, coord2D &klo, coord2D &kup ) const
{
    locate_value( v.x,X, klo.x ,kup.x);
    locate_value( v.y,Y, klo.y, kup.y);
}
