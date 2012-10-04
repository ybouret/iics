#include "cell.hpp"

void Segmenter:: dispatch_vertical_junctions( const mpi &MPI, Cell &cell )
{
    if( MPI.IsFinal )
    {
        //----------------------------------------------------------------------
        // encoding extraneous vertical junctions
        //----------------------------------------------------------------------
        vector<JPack> &jcom = cell.jcom;
        jcom.free();
        assert(0==jcom.size());
        for( unit_t i=X.lower;i<=X.upper;++i)
        {
            const  Segment &seg = Vert(i);
            const  Junction *J  = seg.tail;
            while( J && J->vertex.y >= cell.Y[cell.upper.y] )
            {
                const JPack jpack(i,J);
                jcom.push_back(jpack);
                J = J->prev;
            }
        }
        const size_t count = jcom.size();
        if(count>0)
        {
            fprintf( stderr, "Need to send %lu >= %g\n",count,cell.Y[cell.upper.y] );
            for(size_t i=1;i<=count;++i)
            {
                fprintf( stderr , "\t@%ld: bubble #%u: y=%g\n", jcom[i].i, jcom[i].b, jcom[i].y);
            }
        }
    }
    
    
    
    
}

