#include "cell.hpp"

void Cell:: save_inter( const string &filename ) const
{
    ios::ocstream fp(filename,false);
    for( const Intersection *I = inter.head; I; I=I->next )
    {
        fp("%.15g %.15g\n", I->vertex.x, I->vertex.y );
    }
}


void Cell:: save_inside( const string &filename ) const
{
    vector<V2D> pts;
    collect_inside(pts);
    ios::ocstream fp(filename,false);
    for( size_t i=pts.size();i>0;--i)
    {
        fp("%.15g %.15g\n", pts[i].x, pts[i].y);
    }
}

void Cell:: save_grid( const string &filename ) const
{
    ios::ocstream fp( filename, false );
    for( unit_t j= lower.y; j<= upper.y;++j)
    {
        if( (j&1) )
        {
            for( unit_t i= lower.x; i <= upper.x; ++i )
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
        else
        {
            for( unit_t i= upper.x; i >= lower.x; --i )
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
    }
    
    for( unit_t i=lower.x; i <= upper.x; ++i )
    {
        if( (i&1) )
        {
            for( unit_t j=lower.y; j<=upper.y;++j)
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
        else 
        {
            for( unit_t j=upper.y; j>=lower.y;--j)
            {
                fp("%g %g\n", X[i], Y[j] );
            }
        }
    }
    
}
