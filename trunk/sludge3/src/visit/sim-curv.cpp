#include "simulation.hpp"

visit_handle Simulation:: get_curve( const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if( !parallel )
    {
        if( name == "junctions")
        {
            if(VisIt_CurveData_alloc(&h) != VISIT_ERROR)
            {
                visit_handle hxc, hyc;
                VisIt_VariableData_alloc(&hxc);
                VisIt_VariableData_alloc(&hyc);
                
                const size_t nj = junctions.count_all();
                vector<Real> jx(nj,0);
                vector<Real> jy(nj,0);
                junctions.to_curve(jx, jy);
                
                VisIt_VariableData_setDataD(hxc,VISIT_OWNER_COPY,1,nj,&jx[1]);
                VisIt_VariableData_setDataD(hyc,VISIT_OWNER_COPY,1,nj,&jy[1]);
                VisIt_CurveData_setCoordsXY(h, hxc, hyc);
                
                return h;
            }
        }
    }
    
    for( const Bubble *b = bubbles.head;b;b=b->next)
    {
        const string bn = vformat("bubble%u", unsigned(b->UID));
        if( bn == name )
        {
            if(VisIt_CurveData_alloc(&h) != VISIT_ERROR)
            {
                visit_handle hxc, hyc;
                VisIt_VariableData_alloc(&hxc);
                VisIt_VariableData_alloc(&hyc);
                
                const size_t np = b->size+1;
                vector<Real> bx(np,0);
                vector<Real> by(np,0);
                const Tracer *tr = b->root;
                for(size_t i=1;i<=np;++i,tr=tr->next)
                {
                    bx[i] = tr->pos.x;
                    by[i] = tr->pos.y;
                }
                VisIt_VariableData_setDataD(hxc,VISIT_OWNER_COPY,1,np,&bx[1]);
                VisIt_VariableData_setDataD(hyc,VISIT_OWNER_COPY,1,np,&by[1]);
                VisIt_CurveData_setCoordsXY(h, hxc, hyc);
            }
        }
    }
    return h;
    
}
