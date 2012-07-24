#include "segmenter.hpp"

void Segmenter:: assign_markers()
{
    for( unit_t j=Y.lower;j<=Y.upper;++j)
    {
        Segment::List &Sj = horizontal[j];
        fprintf( stderr, "seg@y= %8.3f: %3lu:", Y[j], Sj.size);
        for( const Segment *s = Sj.head; s; s=s->next )
        {
            fprintf( stderr, " (%g,%g)", s->handle->vertex.x, s->handle->vertex.y);
        }
        fprintf( stderr, "\n");
        
        //unit_t i = X.lower;

    }
}
