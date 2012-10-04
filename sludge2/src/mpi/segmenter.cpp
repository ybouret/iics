#include "cell.hpp"

void Segmenter:: dispatch_vertical_junctions( const mpi &MPI, const Cell &cell )
{
    if( MPI.IsFinal )
    {
        for( unit_t i=X.lower;i<=X.upper;++i)
        {
            const  Segment &seg = Vert(i);
            size_t count        = 0;
            const  Junction *J  = seg.tail;
            while( J && J->vertex.y >= cell.Y[cell.upper.y] )
            {
                ++count;
                J = J->prev;
            }
            if(count>0)
            {
                fprintf( stderr, "@i=%ld, Need to send %lu\n",i,count);
            }
        }
        
    }
    
    
}

