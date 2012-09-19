#include "segmenter.hpp"

void Segmenter:: process( const Bubbles &bubbles )
{
    assert(hseg!=NULL);
    assert(vseg!=NULL);
    for( size_t i=segcount;i>0;--i) segments[i]->empty();
    for( const Bubble *bubble = bubbles.first(); bubble; bubble=bubble->next)
    {
        process_bubble(bubble);
    }
}

void Segmenter:: process_bubble( const Bubble *bubble )
{
    assert(bubble);
    
    //--------------------------------------------------------------------------
    // outer loop on bubble spots
    //--------------------------------------------------------------------------
    for( const Spot *spot=bubble->spots.head;spot;spot=spot->next)
    {
        process_spot(spot);
    }
    
}

void Segmenter:: process_spot( const Spot *spot )
{
    const Tracer *tracer = spot->handle;
    const Vertex  vertex = tracer->vertex;
    assert( vertex.x >= X[X.lower]);
    assert( vertex.x <  X[X.upper]);
    assert( vertex.y >= Y[Y.lower]);
    assert( vertex.y <  Y[Y.upper]);
    
    //--------------------------------------------------------------------------
    // locate the tracer
    //--------------------------------------------------------------------------
    locate_vertex(vertex, spot->klo, spot->kup);
    
    //--------------------------------------------------------------------------
    // compute the absolute next coordinate
    //--------------------------------------------------------------------------
    const Vertex target = vertex + tracer->edge;
    
    //--------------------------------------------------------------------------
    // find junctions
    //--------------------------------------------------------------------------
    compute_junctions(spot, vertex, target);
}

void Segmenter:: compute_junctions( const Spot *spot, const Vertex &self, const Vertex &other )
{
    
}
