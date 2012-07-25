#include "segmenter.hpp"

void Segmenter:: assign_markers()
{
    
    for( unit_t j=Y.lower;j<=Y.upper;++j)
    {
        Segment::List &Sj = horizontal[j];
        if(false)
        {
            fprintf( stderr, "seg@y= %8.3f: %3lu:", Y[j], Sj.size);
            for( const Segment *s = Sj.head; s; s=s->next )
            {
                fprintf( stderr, " (%g,%g)", s->handle->vertex.x, s->handle->vertex.y);
            }
            fprintf( stderr, "\n");
        }
        
        //----------------------------------------------------------------------
        // initialize lookup
        //----------------------------------------------------------------------
        unit_t   i      = X.lower;
        bool     inside = false;    //! TODO: check if bubble has touching points ?
        Segment *s      = Sj.head;
        while( s!=0 )
        {
            //------------------------------------------------------------------
            //-- determine how many times we cross the bubble
            //------------------------------------------------------------------
            size_t num_cross = 1;
            while( (s->next !=0)  && (s->next->handle->lo == s->handle->lo) )
            {
                ++num_cross;
                s=s->next;
            }
            
            //------------------------------------------------------------------
            //-- take action between current i and s->lo
            //------------------------------------------------------------------
            if( inside )
            { 
                const Junction *J = s->handle;
                assert(J->bubble!=NULL);
                assert(J->bubble->id>0);
                while( i <= s->handle->lo )
                {
                    const Coord c(i,j);
                    J->bubble->markers.append()->coord = c;
                    ++i;
                }
            }
            
            //------------------------------------------------------------------
            //-- update status
            //------------------------------------------------------------------
            if( (num_cross&1) )
                inside = !inside;
            
            i=s->handle->up;
            s=s->next;
        }
    }
}
