#include "cell.hpp"

void Segmenter:: dispatch_vertical_junctions( const mpi &MPI, Cell &cell )
{
    static const int tag = 0xD15;
    const int  source = MPI.CommWorldLast;
    const int  target = 0;
    const int  rank   = MPI.CommWorldRank;
    const bool is_source = rank == source;
    const bool is_target = rank == target;
    //const bool is_active = is_source || is_target;
    
    
    MPI_Request      source_request;
    MPI_Request      target_request;
    unsigned long    count = 0;
    vector<JPack>   &jcom  = cell.jcom;
    if( is_source )
    {
        //----------------------------------------------------------------------
        // encoding extraneous vertical junctions
        //----------------------------------------------------------------------
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
        
        //----------------------------------------------------------------------
        // sending #count to target
        //----------------------------------------------------------------------
        count = jcom.size();
        fprintf( stderr, "\t@source: Need to send %lu >= %g\n",count,cell.Y[cell.upper.y] );
        for(size_t i=1;i<=count;++i)
        {
            fprintf( stderr , "\t\t-->@%ld: bubble #%u: y=%g\n", jcom[i].i, jcom[i].b, jcom[i].y);
        }
        MPI.Isend(&count, 1, MPI_UNSIGNED_LONG, target, tag, MPI_COMM_WORLD, source_request);
    }
    
    
    if( is_target )
    {
        //----------------------------------------------------------------------
        // recv #count
        //----------------------------------------------------------------------
        MPI.Irecv(&count, 1, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, target_request);
    }
    
    if( is_source )
    {
        //----------------------------------------------------------------------
        // wait for send of #count
        //----------------------------------------------------------------------
        MPI_Status status;
        MPI.Wait(source_request, status);
        
        //----------------------------------------------------------------------
        // send data to target
        //----------------------------------------------------------------------
        MPI.Isend( jcom(), count * sizeof(JPack), MPI_BYTE, target, tag, MPI_COMM_WORLD, source_request);
    }
    
    if( is_target)
    {
        //----------------------------------------------------------------------
        // wait for recv of #count
        //----------------------------------------------------------------------
        MPI_Status status;
        MPI.Wait(target_request, status);
        
        //----------------------------------------------------------------------
        // ok, now recv data from source
        //----------------------------------------------------------------------
        fprintf( stderr, "\t@target: Need to recv %lu\n",count );
        const JPack invalid_jpack;
        jcom.make( count, invalid_jpack);
        MPI.Irecv( jcom(), count * sizeof(JPack), MPI_BYTE, source, tag, MPI_COMM_WORLD, target_request);
    }
    
    if( is_source )
    {
        MPI_Status status;
        MPI.Wait(source_request, status);
        // source is done
    }
    
    if( is_target )
    {
        MPI_Status status;
        MPI.Wait(target_request, status);
        for( size_t k=1;k<=count;++k)
        {
            const JPack &jpack = jcom[k];
            fprintf( stderr , "\t\t<--@%ld: bubble #%u: y=%g\n", jpack.i, jpack.b, jpack.y);
        }
    }
    
}

