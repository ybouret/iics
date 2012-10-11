#include "../simulation.hpp"


visit_handle Simulation:: get_variable( int domain, const string &name ) const
{

    visit_handle h = VISIT_INVALID_HANDLE;
    
    if( name == "P" )
    {
        const int nComponents= 1;
        const int nTuples    = P.items;
        //MPI.Printf0( stderr, "Sending P: %dx%d\n", nComponents, nTuples);
        assert(P.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, P.entry);
        }
        return h;
    }
    
    if( name == "B" )
    {
        const int nComponents= 1;
        const int nTuples    = B.items;
        //MPI.Printf0( stderr, "Sending B: %dx%d\n", nComponents, nTuples);
        assert(B.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, B.entry);
        }
        return h;
    }
    
    if( name == "Bulk" )
    {
        const int nComponents= 1;
        const int nTuples    = B.items;
        //MPI.Printf0( stderr, "Sending B: %dx%d\n", nComponents, nTuples);
        assert(Bulk.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, Bulk.entry);
        }
        return h;
    }

    
    if( name == "U" )
    {
        const int nComponents= 2;
        const int nTuples    = U.items;
        //MPI.Printf0( stderr, "Sending U: %dx%d\n", nComponents, nTuples);
        assert(U.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real*)(U.entry));
        }
        return h;
    }
    
    if( name == "gradP" )
    {
        const int nComponents= 2;
        const int nTuples    = gradP.items;
        //MPI.Printf0( stderr, "Sending U: %dx%d\n", nComponents, nTuples);
        assert(gradP.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real*)(gradP.entry));
        }
        return h;
    }

    if( name == "Penter" )
    {
        const int nComponents= 2;
        const int nTuples    = U.items;
        assert(Penter.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real*)(Penter.entry));
        }
        return h;
    }
    
    if( name == "Pleave" )
    {
        const int nComponents= 2;
        const int nTuples    = U.items;
        assert(Pleave.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real*)(Pleave.entry));
        }
        return h;
    }


    return h;
}

#include "yocto/string/conv.hpp"

visit_handle Simulation:: get_curve( const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    
    if( strncmp(name.c_str(),"bubble",6) == 0 )
    {
        const size_t  id = strconv::to_size( name.c_str()+6, "bubble id");
        assert(id>0);
        assert(id<=bubbles.count());
        const Bubble *b  = bubbles.first();
        for(size_t iter=1;iter<id;++iter) b=b->next;
        assert(b);
        if( VisIt_CurveData_alloc(&h) == VISIT_OKAY )
        {
            // copy bubble spots coordinates
            const size_t bn = b->size+1;
            vector<Real> bx(bn,0);
            vector<Real> by(bn,0);
            const Tracer *p = b->root;
            
            
            for( size_t i=1; i <= bn; ++i,p=p->next)
            {
                bx[i] = p->vertex.x;
                by[i] = p->vertex.y;
            }
            
            
            // make a curve
            visit_handle hcx,hcy;
            VisIt_VariableData_alloc( &hcx );
            VisIt_VariableData_alloc( &hcy );
            VisIt_VariableData_setDataD(hcx, VISIT_OWNER_COPY, 1, bn, bx());
            VisIt_VariableData_setDataD(hcy, VISIT_OWNER_COPY, 1, bn, by());
            VisIt_CurveData_setCoordsXY(h, hcx, hcy);
            
        }
    }
    
    if( name == "junctions" )
    {
        if( VisIt_CurveData_alloc(&h) == VISIT_OKAY )
        {
            // copy juntions coordinates
            const size_t nj = segmenter.num_junctions();
            vector<Real> jx(nj,0);
            vector<Real> jy(nj,0);
            const Segments &segments = segmenter();
            for( size_t i=segments.size(),j=1;i>0;--i)
            {
                const Segment &seg = *segments[i];
                for( const Junction *J = seg.head;J;J=J->next)
                {
                    assert(j<=nj);
                    jx[j] = J->vertex.x;
                    jy[j] = J->vertex.y;
                    ++j;
                }
            }
            
            // make a curve
            visit_handle hcx,hcy;
            VisIt_VariableData_alloc( &hcx );
            VisIt_VariableData_alloc( &hcy );
            VisIt_VariableData_setDataD(hcx, VISIT_OWNER_COPY, 1, nj, jx());
            VisIt_VariableData_setDataD(hcy, VISIT_OWNER_COPY, 1, nj, jy());
            VisIt_CurveData_setCoordsXY(h, hcx, hcy);

        }

    }
    
    return h;
}
