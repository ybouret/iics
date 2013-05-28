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
	const size_t np = max_of<size_t>( 3, numeric<Real>::two_pi * radius / b->lambda );
	b->auto_delete();
	for(size_t i=0; i < np; ++i )
	{
		const Real theta = (numeric<Real>::two_pi * i)/np;
		Tracer *   tr    = new Tracer();
		tr->pos.x = radius * Cos( theta ) + center.x;
		tr->pos.y = radius * Sin( theta ) + center.y;
		b->push_back(tr);
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
    b->auto_delete();

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
        Tracer *   tr    = new Tracer();
		tr->pos.x = R.x * Cos( t ) + C.x;
		tr->pos.y = R.y * Sin( t ) + C.y;
		b->push_back(tr);
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
}

