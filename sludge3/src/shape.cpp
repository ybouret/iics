#include "shape.hpp"
#include "yocto/code/utils.hpp"

void Shape::Circle( Bubble *b, const Vertex center, Real radius ) throw()
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
