#include "bubble.hpp"
#include "yocto/ios/ocstream.hpp"

void Bubble:: save_dat(ios::ostream &fp) const
{
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

void Bubble:: save_dat( const string &fn ) const
{
    ios::ocstream fp(fn,false);
    save_dat(fp);
}