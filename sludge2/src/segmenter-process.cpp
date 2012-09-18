#include "segmenter.hpp"

void Segmenter:: process( const Bubbles &bubbles )
{
    assert(hseg!=NULL);
    assert(vseg!=NULL);
    for( size_t i=segcount;i>0;--i) segments[i]->empty();
    for( const Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next)
    {
        process1(bubble);
    }
}

void Segmenter:: process1( const Bubble *bubble )
{
    assert(bubble);
    
    //--------------------------------------------------------------------------
    // outer loop on bubble spots
    //--------------------------------------------------------------------------
    for( const Spot *spot=bubble->spots.head;spot;spot=spot->next)
    {
        const Tracer *tracer = spot->handle;
        const Vertex  vertex = tracer->vertex;
        assert( vertex.x >= X[X.lower]);
        assert( vertex.x <  X[X.upper]);
        assert( vertex.y >= Y[Y.lower]);
        assert( vertex.y <  Y[Y.upper]);
        
        //----------------------------------------------------------------------
        // locate the tracer
        //----------------------------------------------------------------------
        locate_vertex(vertex, spot->klo, spot->kup);
        
        
        
        
    }
    
}