#include "shape.hpp"
#include "yocto/code/utils.hpp"

////////////////////////////////////////////////////////////////////////////////
//
// CIRCLE
//
////////////////////////////////////////////////////////////////////////////////
void Shape::Circle( Bubble *b, const Vertex center, Real radius )
{
	assert(b);
	assert(radius>0);
    
	b->clear();
    
    const size_t np = max_of<size_t>( 3, numeric<Real>::two_pi * radius / b->lambda );
    for(size_t i=0; i < np; ++i )
	{
		const Real theta = (numeric<Real>::two_pi * i)/np;
		Tracer *   tr    = b->append();
		tr->pos.x = radius * Cos( theta ) + center.x;
		tr->pos.y = radius * Sin( theta ) + center.y;
	}
    
	b->init_contour();
}


////////////////////////////////////////////////////////////////////////////////
//
// Ellipse
//
////////////////////////////////////////////////////////////////////////////////
#include "yocto/sequence/vector.hpp"

namespace
{
    class MakeEllipse
    {
    public:
        const size_t np;
        vector<Real> theta;
        const Vertex R;
        
        MakeEllipse( size_t n , const Vertex &r) :
        np(n+1),
        theta(np,0),
        R(r)
        {
            
        }
        
        ~MakeEllipse() throw()
        {
        }
        
        Real build( Real mu ) throw()
        {
            theta[1] = 0;
            for(size_t i=2;i<=np;++i)
            {
                const Real t = theta[i-1];
                const Real sa = R.x * Sin(t);
                const Real cb = R.y * Cos(t);
                theta[i] = theta[i-1] + mu / Sqrt(sa*sa+cb*cb);
            }
            return theta[np] - numeric<Real>::two_pi;
        }
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(MakeEllipse);
    };
}

#include "yocto/math/fcn/zfind.hpp"

void Shape:: Ellipse( Bubble *b, const Vertex C, const Vertex R)
{
    assert(b);
    assert(R.x>0);
    assert(R.y>0);
    
    b->clear();

    //==========================================================================
    // approx perimeter
    //==========================================================================
    const Real   A = R.x;
    const Real   B = R.y;
    const Real   AmB = A-B;
    const Real   L = numeric<Real>::pi * Sqrt( 2.0 * (A*A+B*B) - AmB*AmB/2.0);
    
    //==========================================================================
    // build approx theta
    //==========================================================================
    MakeEllipse ell(max_of<size_t>( 3, L/b->lambda),R);
    numeric<Real>::function zell( &ell, &MakeEllipse::build);
    assert(zell(0)<=0);
    Real mu = b->lambda;
    while( zell(mu) <=0 ) mu *= 1.1;
    zfind<Real> solve(1e-4);
    mu = solve(zell,0,mu);
    
    //==========================================================================
    // build points
    //==========================================================================
    for( size_t i=1; i<ell.np;++i)
    {
        const Real t = ell.theta[i];
        Tracer *   tr    = b->append();
		tr->pos.x = R.x * Cos( t ) + C.x;
		tr->pos.y = R.y * Sin( t ) + C.y;
    }
    
    b->init_contour();
    
    
}


////////////////////////////////////////////////////////////////////////////////
//
// rotate
//
////////////////////////////////////////////////////////////////////////////////
void Shape:: Rotate(Bubble *b, const Real alpha)
{
    assert(b);
    assert(b->size>=3);
    Tracer *tr = b->root;
    const Real ca = Cos(alpha);
    const Real sa = Sin(alpha);
    
    for(size_t i=b->size;i>0;--i,tr=tr->next)
    {
        const Vertex ray(b->G,tr->pos);
        Vertex r( ray.x * ca - ray.y * sa, ray.x *sa + ray.y *ca);
        tr->pos = b->G + r;
    }
    
    b->init_contour();
    
    
}

////////////////////////////////////////////////////////////////////////////////
//
// Blob
//
////////////////////////////////////////////////////////////////////////////////
void Shape:: Blob( Bubble *b, const Vertex C, const Real radius,  Real rho,  Real w)
{
    assert(b);
    assert(radius>0);
    b->clear();
    
    rho = clamp<Real>(0,rho,0.95);
    w   = clamp<Real>(0,w,1);
    const Real   fac = radius/(1.0 + rho);
    const Real   L   = 20 * radius;
    const size_t np = max_of<size_t>(3,L/b->lambda);
    
    const Real A = rho * (1-w);
    const Real B = rho * w;
    
    for( size_t i=0; i < np; ++i )
    {
        const Real theta = (numeric<Real>::two_pi * i)/np;
		Tracer *   tr    = b->append();
        const Real rp    = fac * ( 1.0 + A * Cos( 2*theta) + B * Cos(3*theta) );
		tr->pos.x = rp * Cos( theta ) + C.x;
		tr->pos.y = rp * Sin( theta ) + C.y;
    }
    b->init_contour();
    
}

////////////////////////////////////////////////////////////////////////////////
//
// Grow
//
////////////////////////////////////////////////////////////////////////////////

void Shape:: Grow( Bubble *b, const Real factor)
{
    assert(b);
    assert(b->size>=3);
    Tracer *tr = b->root;
    for(size_t i=b->size;i>0;--i,tr=tr->next)
    {
        const Vertex r(b->G,tr->pos);
        tr->pos = b->G + factor * r;
    }
    b->init_contour();
    
}

////////////////////////////////////////////////////////////////////////////////
//
// Square
//
////////////////////////////////////////////////////////////////////////////////
void Shape:: Square( Bubble *b, const Vertex C, Real a)
{
    assert(b);
    assert(a>0);
    b->clear();
    
    const size_t nextra = size_t(ceil(a/b->lambda));

    const Real h = a/2;
    const Real step = a/(nextra+1);
    
    // bottom|left -> bottom|right
    Vertex org(-h,-h);
    org += C;
    b->append(org);
    for(size_t i=1; i <= nextra; ++i)
    {
        const Vertex dv(step*i,0);
        b->append( org + dv );
    }
    
    // bottom|right -> top/right
    org = Vertex(h,-h);
    org += C;
    b->append(org);
    for(size_t i=1; i <= nextra; ++i)
    {
        const Vertex dv(0,step*i);
        b->append( org + dv );
    }
    
    // top|right -> top|left
    org = Vertex(h,h);
    org += C;
    b->append(org);
    for(size_t i=1; i <= nextra; ++i)
    {
        const Vertex dv(step*i,0);
        b->append( org - dv );
    }

    // top|left -> bottom|left
    org = Vertex(-h,h);
    org += C;
    b->append(org);
    for(size_t i=1; i <= nextra; ++i)
    {
        const Vertex dv(0,step*i);
        b->append( org - dv );
    }
    
    b->init_contour();
    
}

void Shape:: Move( Bubble *b, const Vertex v)
{
    assert(b);
    Tracer *tr = b->root;
    for(size_t i=b->size;i>0;--i,tr=tr->next)
    {
        tr->pos += v;
    }
    b->init_contour();
}


