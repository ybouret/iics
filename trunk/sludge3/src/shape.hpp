#ifndef SLUDGE_SHAPE_INCLUDED
#define SLUDGE_SHAPE_INCLUDED 1

#include "bubble.hpp"

struct Shape
{
	static void Circle( Bubble *b, const Vertex center, Real radius );
    static void Ellipse( Bubble *b, const Vertex C, const Vertex R);
    static void Rotate( Bubble *b, const Real alpha);
    static void Blob( Bubble *b, const Vertex C, Real radius, Real rho, Real w=0);
    static void Square( Bubble *b, const Vertex C, Real a);
    static void Grow( Bubble *b, const Real factor);
    static void Move( Bubble *b, const Vertex v);
    
};

#endif
