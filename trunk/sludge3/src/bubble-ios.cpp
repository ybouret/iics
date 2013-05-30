#include "bubble.hpp"
#include "yocto/ios/ocstream.hpp"


void Bubble:: save_dat( const string &fn ) const
{
    ios::ocstream fp(fn,false);
    if(size)
    {
        const Tracer *tr = root;
        for(size_t i=size;i>0;--i,tr=tr->next)
        {
            fp("%g %g\n", tr->pos.x, tr->pos.y);
        }
        fp("%g %g\n", tr->pos.x, tr->pos.y);
    }
    
}

void Bubble:: save_t( const string &fn ) const
{
    ios::ocstream fp(fn,false);
    if(size)
    {
        const Tracer *tr = root;
        for(size_t i=size;i>0;--i,tr=tr->next)
        {
            const Vertex &p = tr->pos;
            const Vertex  q = p + (lambda * tr->t);
            
            fp("%g %g\n",p.x,p.y);
            fp("%g %g\n",q.x,q.y);
            fp("\n");
        }
    }
}

void Bubble:: save_n( const string &fn ) const
{
    ios::ocstream fp(fn,false);
    if(size)
    {
        const Tracer *tr = root;
        for(size_t i=size;i>0;--i,tr=tr->next)
        {
            const Vertex &p = tr->pos;
            const Vertex  q = p + (tr->C * lambda * tr->n);
                        
            fp("%g %g\n",p.x,p.y);
            fp("%g %g\n",q.x,q.y);
            fp("\n");
        }
    }
}


void Bubble:: save_all(const string &pfx) const
{
    {
        const string fn = pfx + ".dat";
        save_dat(fn);
    }
    
    {
        const string fn = pfx + "_t.dat";
        save_t(fn);
    }
    
    {
        const string fn = pfx  + "_n.dat";
        save_n(fn);
    }
    
}

