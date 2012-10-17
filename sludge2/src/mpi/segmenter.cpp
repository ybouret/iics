#include "cell.hpp"

void Segmenter:: dispatch_vertical_junctions( const mpi &MPI, Cell &cell )
{
    static const int tag = 0xD15;
    const int  source    = MPI.CommWorldLast;
    const int  target    = 0;
    const int  rank      = MPI.CommWorldRank;
    const bool is_source = (rank == source);
    const bool is_target = (rank == target);
    
    
    MPI_Request      source_request;
    MPI_Request      target_request;
    unsigned long    count = 0;
    vector<JPack>   &jsend = cell.jsend;
    vector<JPack>   &jrecv = cell.jrecv;
    if( is_source )
    {
        //----------------------------------------------------------------------
        // encoding extraneous vertical junctions
        //----------------------------------------------------------------------
        const Real ylo = cell.Y[cell.upper.y-1];
        const Real yhi = cell.Y[cell.upper.y+1];
        jsend.free();
        assert(0==jsend.size());
        for( unit_t i=X.lower;i<=X.upper;++i)
        {
            const  Segment &seg = Vert(i);
            const  Junction *J  = seg.tail;
            while( J && J->vertex.y >=  ylo && J->vertex.y < yhi)
            {
                const JPack jpack(i,J);
                jsend.push_back(jpack);
                J = J->prev;
            }
        }
        
        //----------------------------------------------------------------------
        // sending #count to target
        //----------------------------------------------------------------------
        count = jsend.size();
#if 0
        fprintf( stderr, "\t@source: Need to send %lu\n",count);
        for(size_t i=1;i<=count;++i)
        {
            const JPack &jpack = jsend[i];
            fprintf( stderr , "\t\t-->@%ld: bubble #%u: x=%g, y=%g, c=%g\n", jpack.i, jpack.b, X[jpack.i], jpack.y, jpack.c);
        }
#endif
        MPI.Isend(&count, 1, MPI_UNSIGNED_LONG, target, tag, MPI_COMM_WORLD, source_request);
    }
    
    
    if( is_target )
    {
        remove_vertical_junctions_below(cell.pbc.lo);
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
#if !defined(NDEBUG)
        if( count > 0 )
        {
            assert( jsend() == &jsend[1] );
        }
#endif
        MPI.Isend( jsend(), count * sizeof(JPack), MPI_BYTE, target, tag, MPI_COMM_WORLD, source_request);
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
        //fprintf( stderr, "\t@target: Need to recv %lu\n",count );
        const JPack invalid_jpack;
        jrecv.make( count, invalid_jpack);
#if !defined(NDEBUG)
        if( count > 0 )
        {
            assert(jrecv() == &jrecv[1]);
        }
#endif
        MPI.Recv( jrecv(), count * sizeof(JPack), MPI_BYTE, source, tag, MPI_COMM_WORLD, status);
        
        //----------------------------------------------------------------------
        // unpacking
        //----------------------------------------------------------------------
        for( size_t k=1;k<=count;++k)
        {
            const JPack &jpack = jrecv[k];
            //fprintf( stderr , "\t\t<--@%ld: bubble #%u: x=%g, y=%g, c=%g\n", jpack.i, jpack.b, X[jpack.i], jpack.y-cell.pbc.L,jpack.c);
            Junction *J  = Vert(jpack.i).append();
            J->kind      = Junction::Vert;
            J->vertex.x  = X[jpack.i];
            J->vertex.y  = jpack.y - cell.pbc.L;
            J->curvature = jpack.c;
            J->n         = jpack.n;
            J->t         = jpack.t;
            J->pressure  = jpack.p;
            J->gradP_t   = jpack.g;
            J->bubble    = cell.bubbles.first();
            assert(jpack.b>0);
            for(size_t ii=jpack.b;ii>1;--ii)
            {
                J->bubble = J->bubble->next;
                assert(J->bubble);
            }
            locate_value( J->vertex.y, Y, J->klo, J->khi);
        }
        SortVert();
        // target is done
    }
    
    if( is_source )
    {
        MPI_Status status;
        MPI.Wait(source_request, status);
        // source is done
    }
    
}

