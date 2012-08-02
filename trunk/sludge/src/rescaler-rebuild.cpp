#include "rescaler.hpp"
#include "yocto/ios/ocstream.hpp"

//! s[i] <= sx < s[i+1]
static inline
size_t __locate_abscissa( const Real sx, const array<Real> &s )
{
    const size_t n = s.size();
    assert( sx <  s[n] );
    assert( sx >= s[1] );
    size_t ilo = 1;
    size_t iup = n;
    while( (iup-ilo) > 1 )
    {
        const size_t mid = (ilo+iup)>>1; assert(mid>0);
        if( sx < s[mid])
        {
            iup = mid;
        }
        else
        {
            ilo = mid;
        }
    }
    assert(ilo>0);
    assert(iup<=n);
    return ilo;
}

static inline
Vertex __interpv( Real sx, const Real sa[], const Real xa[], const Real ya[], size_t n )
{
    Vertex v(0,0);
    for( size_t i=0; i <n; ++i )
    {
        const Real si = sa[i];
        Real  num=1;
        Real  den=1;
        for( size_t j=0; j<i;++j )
        {
            const Real sj = sa[j];
            num  *= sx - sj;
            den  *= si - sj;
        }
        for( size_t j=i+1; j<n;++j )
        {
            const Real sj = sa[j];
            num  *= sx - sj;
            den  *= si - sj;
        }
        v.x += (num*xa[i])/den;
        v.y += (num*ya[i])/den;
    }
    return v;
}

void Rescaler:: rebuild( Bubble &bubble )
{
    //static int fid = 0;
    
    const size_t n  = bubble.size;
    const size_t n1 = n+1;
    assert( s.size()     == n1 );
    assert( ax.size()    == n1 );
    assert( ay.size()    == n1 );
    assert(period>0);
    assert(bubble.area>0);
    assert(a_list.size >= 3);
    
#if 0
    {
        ios::ocstream fp( vformat("org%d.dat",fid), false);
        for( size_t i=1; i <= n1; ++i )
        {
            fp("%g %g %g\n", s[i], ax[i], ay[i]);
        }
        fp("%g %g %g\n", period, ax[1], ay[1]);
    }
#endif
    
    Real sa[4] = { 0 };
    Real xa[4] = { 0 };
    Real ya[4] = { 0 };
    
    //--------------------------------------------------------------------------
    // rebuild the bubble
    //--------------------------------------------------------------------------
    //std::cerr << "rebuilding " << a_list.size << " points" << std::endl;
    bubble.empty();
    for( const abscissa *a = a_list.head; a; a=a->next )
    {
        const Real s_i = a->s;
        assert(s_i>=0);
        assert(s_i<period);
        //-- locate the abscissa
        const size_t j1  = __locate_abscissa(s_i, s);
        const size_t j2  = j1+1;
        
        //-- fill the middle points
        sa[1] = s[j1]; xa[1] = ax[j1]; ya[1] = ay[j1];
        sa[2] = s[j2]; xa[2] = ax[j2]; ya[2] = ay[j2];
        
        //-- fill the left point
        if( j1 <= 1 )
        {
            sa[0] = sa[1] - (s[n1] - s[n]); xa[0] = ax[n]; ya[0] = ay[n];
        }
        else
        {
            const size_t j0 = j1-1; assert(j0>0);
            sa[0] = s[j0]; xa[0] = ax[j0]; ya[0] = ay[j0];
        }
        
        //-- fill the right point
        if( j2 >= n1 )
        {
            sa[3] = s[n1] + (s[2]-s[1]);
            xa[3] = ax[1];
            ya[3] = ay[1];
        }
        else
        {
            const size_t j3 = j2+1;
            sa[3] = s[j3]; xa[3] = ax[j3]; ya[3] = ay[j3];
        }
#if 1
        const Vertex v = __interpv(s_i, sa, xa, ya, 4);
#else
        const Real    fac = (s_i - sa[1])/(sa[2]-sa[1]);
        const Vertex  v( xa[1] + fac * ( xa[2] - xa[1]), ya[1] + fac * ( ya[2] - ya[1]) );
#endif
        bubble.append()->vertex = v;
    }
    assert(bubble.size==a_list.size);
    
#if 0
    {
        ios::ocstream fp( vformat("ref%d.dat",fid), false);
        const abscissa *a = a_list.head;
        const Tracer   *p = bubble.root;
        while(a)
        {
            fp("%g %g %g\n", a->s, p->vertex.x, p->vertex.y);
            a=a->next;
            p=p->next;
        }
        fp("%g %g %g\n", period, p->vertex.x, p->vertex.y);
        
        ++fid;
    }
#endif
    
    //--------------------------------------------------------------------------
    // rebuild its metrics
    //--------------------------------------------------------------------------
    build_metrics(bubble,RescaleWithConstantPressure);
    
}
