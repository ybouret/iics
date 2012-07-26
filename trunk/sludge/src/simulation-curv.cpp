#include "simulation.hpp"

visit_handle Simulation:: get_curve( const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if( name == "bubble" )
    {
        const Bubble *b = bubbles.first();
        assert(b);
        if( VisIt_CurveData_alloc(&h) == VISIT_OKAY )
        {
            // copy bubble spots coordinates
            const size_t bn = b->size+1;
            vector<Real> bx(bn,0);
            vector<Real> by(bn,0);
            const Tracer *p = b->root;
            
            //ios::ocstream fp("bubble.dat",false);
            
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
    return h;
}
