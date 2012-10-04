#include "cell.hpp"

void Segmenter:: dispatch_vertical_junctions( const mpi &MPI, Cell &cell )
{
    static const int tag = 0xD15;
    const int source = MPI.CommWorldLast;
    const int target = 0;
    const int rank   = MPI.CommWorldRank;
    const bool is_source = rank == source;
    const bool is_target = rank == target;
    
    
    MPI_Request      request;
    unsigned long    count = 0;
    if( is_source )
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
        count = jcom.size();
        fprintf( stderr, "\t@source: Need to send %lu >= %g\n",count,cell.Y[cell.upper.y] );
        for(size_t i=1;i<=count;++i)
        {
            fprintf( stderr , "\t\t@%ld: bubble #%u: y=%g\n", jcom[i].i, jcom[i].b, jcom[i].y);
        }
        MPI.Isend(&count, 1, MPI_UNSIGNED_LONG, target, tag, MPI_COMM_WORLD, request);
    }
    
    if( is_target )
    {
        MPI.Irecv(&count, 1, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, request);
    }
    
    if( is_source || is_target )
    {
        MPI_Status status;
        MPI.Wait(request, status);
    }
    
    if( is_target)
    {
        fprintf( stderr, "\t@target: Need to recv %lu\n",count );
    }
    
    
}

