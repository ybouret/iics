#include "fish.hpp"
#include "yocto/ios/ocstream.hpp"

void Fish:: save_vtk( const string &filename ) const
{
    const size_t NP = points.size();
    ios::ocstream fp(filename,false);
    fp << "# vtk DataFile Version 3.0\n";
    fp << "Fish\n";
    fp << "ASCII\n";
    fp << "DATASET POLYDATA\n";
    fp("POINTS %u float\n", unsigned(NP));
    for(size_t i=1;i<=NP;++i)
    {
        const vtx_t &r = points[i]->r;
        fp("%g %g %g\n",r.x,r.y,r.z);
    }

    const size_t NT = triangles.size();

    fp("POLYGONS %u %u\n", unsigned(NT), unsigned(4*NT) );
    for(size_t i=1;i<=NT;++i)
    {
        const Triangle &tr = triangles[i];
        fp("3 %u %u %u\n", unsigned(tr.a->i), unsigned(tr.b->i), unsigned(tr.c->i) );
    }

}

inline void __stl( ios::ostream &fp, const Point &p )
{
    fp("   vertex %.6g %.6g %.6g\n", p.r.x, p.r.y,p.r.z);
}

void Fish:: save_stl( const string &filename ) const
{
    ios::ocstream fp(filename,false);

    fp("solid fish\n");

    const size_t NT = triangles.size();
    for(size_t i=1;i<=NT;++i)
    {
        const Triangle &tr = triangles[i];
        fp(" facet normal %.6g %.6g %.6g\n",tr.n.x,tr.n.y,tr.n.z);
        fp("  outer loop\n");
        __stl(fp, *tr.a);
        __stl(fp, *tr.b);
        __stl(fp, *tr.c);
        fp("  end loop\n");
        fp(" end facet\n");
    }
    fp("endsolid\n");
}

