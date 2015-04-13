#include "fish.hpp"
#include "yocto/exception.hpp"

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
CubiX:: ~CubiX() throw() {}

static inline
void LoadCoordinates( sequence<vtx_t> &beta, lua_State *L, const string &name)
{
    lua_settop(L, 0);
    lua_getglobal(L,name.c_str());
    if( !lua_istable(L, -1) )
        throw exception("No Table '%s'", name.c_str());

    const size_t n = lua_rawlen(L, -1);
    if(n<1)
        throw exception("Need At Least 1 point in '%s'", name.c_str());

    beta.free();
    beta.ensure(n+2);

    {
        vtx_t tmp(0,0,0);
        beta.push_back(tmp);
    }
    for(size_t i=1;i<=n;++i)
    {
        lua_rawgeti(L, -1, i);
        if(!lua_isnumber(L, -1))
        {
            throw exception("%s[%u] is not a number", name.c_str(), unsigned(i));
        }
        const double t = lua_tonumber(L, -1);
        lua_pop(L,1);
        vtx_t v(t,t,0);

        if(i<=1)
        {
            v.z=0;
        }
        else
        {
            if(i>=n)
            {
                v.z = 1;
            }
            else
            {
                v.z = double(i-1)/(n-1);
            }
        }
        beta.push_back(v);
    }

    {
        vtx_t tmp(0,0,1);
        beta.push_back(tmp);
    }
    
    
}


CubiX::CubiX( lua_State *L, const string &table_name, const size_t iCoord ) :
CubiXBase(),
solv(0),
zfcn(this, & CubiXBase::Z, 0),
indx(clamp<size_t>(0,iCoord,1))
{
    LoadCoordinates(*this, L, table_name);
    assert(size()>=3);
}

double CubiX:: GetValue( const double z )
{
    assert(size()>=3);
    vtx_t result;
    if(z<=0)
    {
        result = front();
    }
    else
    {
        if(z>=1)
        {
            result = back();
        }
        else
        {
            zfcn.target = z;
            result = Compute(solv(zfcn.call,0,1));
        }
    }
    const double *q = &result.x;
    return q[indx];
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

Profile:: ~Profile() throw() {}


#include "yocto/math/opt/minimize.hpp"

Profile:: Profile( lua_State *L ) :
W(  &width,  & CubiX::GetValue ),
H(  &height, & CubiX::GetValue ),
width( L,"width", 0),
height(L,"height",1),
zarr(NZ,0),
rmax(NZ,0),
arcL(NZ,0),
maxL(0),
maxP(0)
{
    for(size_t i=1;i<=NZ;++i)
    {
        const double zz = double(i-1)/(NZ-1);
        zarr[i]  = zz;
        rmax[i]  = max_of<double>( Fabs( W(zz)), Fabs( H(zz) ) );
    }
    arcL[1] = 0;
    for(size_t i=2;i<=NZ;++i)
    {
        const double dr = rmax[i] - rmax[i-1];
        const double dz = zarr[i] - zarr[i-1];
        arcL[i] = arcL[i-1] + Hypotenuse(dr, dz);
    }
    (double &)maxL = arcL[NZ];
    for(size_t i=2;i<NZ;++i) { arcL[i] /= maxL; }
    arcL[NZ] = 1;

    if(maxL<=0)
        throw exception("Invalid Profile");

    std::cerr << "Max Arc Length=" << maxL << std::endl;

    numeric<double>::function optP(this, &Profile::minusPerimeter);
    triplet<double> ZZ = {0, 0.5, 1.0};
    triplet<double> PP = { optP(ZZ.a), optP(ZZ.b), optP(ZZ.c) };
    minimize<double>(optP, ZZ, PP, 0);

    (double &)maxP = computePerimeter(ZZ.b);
    std::cerr << "Max Perimeter=" << maxP << " @" << ZZ.b << std::endl;
}

double Profile:: getZ( const double ratio )
{
    if(ratio<=0)
    {
        return 0;
    }
    else
    {
        if(ratio>=1)
        {
            return 1;
        }
        else
        {
            size_t jlo = 1;
            size_t jhi = NZ;
            while(jhi-jlo>1)
            {
                const size_t jmid = (jlo+jhi)>>1;
                const double amid = arcL[jmid];
                if(amid<ratio)
                {
                    jlo = jmid;
                }
                else
                {
                    jhi = jmid;
                }
            }
            return zarr[jlo] + (zarr[jhi]-zarr[jlo]) * (ratio-arcL[jlo]) / (arcL[jhi]-arcL[jlo]);
        }
    }
}

double Profile:: computePerimeter( const double z )
{
    const double a = W(z);
    const double b = H(z);
    const double h2 = a*a + b*b;
    return numeric<double>::pi * sqrt(h2+h2);
}


double Profile::minusPerimeter(const double z)
{
    return -computePerimeter(z);
}



////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////

Point::  Point() throw() : i(Index++), r() {}
Point:: ~Point() throw() {}

size_t Point::Index = 0;

////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Slice:: ~Slice() throw() {}

Slice:: Slice(double zz) throw() :
z(zz),
points()
{
}


////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////
Triangle:: ~Triangle() throw() {}


Triangle:: Triangle(const pPoint &A,
                    const pPoint &B,
                    const pPoint &C) throw():
a(A),b(B),c(C),
G( (1.0/3) * (a->r+b->r+c->r) ),
n(),
S(0)
{
    assert(a->i!=b->i);
    assert(a->i!=c->i);
    assert(b->i!=c->i);
    {
        //const vtx_t G = (1.0/3) * (a->r+b->r+c->r);
        const vtx_t AB(a->r,b->r);
        const vtx_t AC(a->r,c->r);
        const vtx_t NN = vtx_t::cross_(AB, AC);
        const vtx_t middle(0,0,0.5);
        const vtx_t Q(middle,G);
        if( (NN*Q) <= 0 )
        {
            b.swap_with(c);
        }
    }

    {
        const vtx_t AB(a->r,b->r);
        const vtx_t AC(a->r,c->r);
        n = vtx_t::cross_(AB, AC);
        S = n.norm();
        n.normalize();
    }
}


Triangle:: Triangle(const Triangle &other) throw() :
a(other.a), b(other.b), c(other.c),
G(other.G),
n(other.n),
S(other.S)
{
}

void Triangle:: recompute() throw()
{
    G = (1.0/3) * (a->r+b->r+c->r);
    const vtx_t AB(a->r,b->r);
    const vtx_t AC(a->r,c->r);
    n = vtx_t::cross_(AB, AC);
    S = n.norm();
    n.normalize();
}


////////////////////////////////////////////////////////////////////////////////
//
// Fish
//
////////////////////////////////////////////////////////////////////////////////
Fish:: Fish( lua_State * L ) : Profile(L)
{
}


Fish:: ~Fish() throw()
{
}

void Fish:: clear() throw()
{
    triangles.release();
    points.release();
    slices.release();
    Point::Index = 0;
}

void Fish:: generateShell( size_t N )
{

    // clean up
    clear();

    // compute the ratio steps
    N = max_of<size_t>(1,N);
    const double delta = 1.0/(N);

    // compute the slices position
    for(size_t i=1;i<=N;++i)
    {
        const double ratio = i*delta;
        const pSlice pS( new Slice( getZ(ratio) ) );
        slices.push_back(pS);
    }

    // compute the number of points per slice M = 2*n+2;
    const size_t n = max_of<size_t>(2,ceil(maxP/delta))-2;
    const size_t M = 2+2*n;
    std::cerr << "\tn=" << n << ", M=" << M << std::endl;
    

    // head point
    pPoint p0( new Point() );
    points.push_back(p0);

    // body points
    for(size_t i=1;i<=N;++i)
    {

        Slice &slice = *slices[i];

        const double w = W(slice.z);
        const double h = H(slice.z);

        //std::cerr << "w=" << width << ", h=" << height << std::endl;
        for(size_t j=0;j<M;++j)
        {
            const double theta = (j*numeric<double>::two_pi)/M;
            pPoint pp( new Point() );

            pp->r.x = w  * cos(theta);
            pp->r.y = h * sin(theta);
            pp->r.z = slice.z;

            slice.points.push_back(pp);
            points.push_back(pp);
        }

    }


    // tail points
    pPoint pN( new Point() );
    points.push_back(pN);
    pN->r.z = 1;

    std::cerr << "#points=" << points.size() << std::endl;

    // generate triangles

    std::cerr << "-- Computing Triangles" << std::endl;
    std::cerr << "\t Head..." << std::endl;

    // head
    {
        const Slice &slice = *slices[1];
        for(size_t i=1;i<=M;++i)
        {
            size_t   ip = i+1;
            if(ip>M) ip = 1;
            const Triangle tr(p0,slice.points[i],slice.points[ip]);
            triangles.push_back(tr);
        }
    }

    std::cerr << "\t Body..." << std::endl;
    // inside
    for(size_t j=1;j<N;++j)
    {
        const array<pPoint> &P0 = slices[j]->points;
        const array<pPoint> &P1 = slices[j+1]->points;

        // loop over quads
        for(size_t i=1;i<=M;++i)
        {
            size_t   ip = i+1;
            if(ip>M) ip = 1;
            const pPoint  &P00 = P0[i];
            const pPoint  &P01 = P0[ip];
            const pPoint  &P10 = P1[i];
            const pPoint  &P11 = P1[ip];

            {
                const Triangle tr(P00,P01,P11);
                triangles.push_back(tr);
            }

            {
                const Triangle tr(P00,P10,P11);
                triangles.push_back(tr);
            }


        }
    }

    std::cerr << "\t Tail..." << std::endl;
    // tail
    {
        const Slice &slice = *slices[N];
        for(size_t i=1;i<=M;++i)
        {
            size_t   ip = i+1;
            if(ip>M) ip = 1;
            const Triangle tr(pN,slice.points[i],slice.points[ip]);
            triangles.push_back(tr);
        }
    }

    std::cerr << "#triangles=" << triangles.size() << std::endl;


}

void Fish:: centerAndRescaleBy( double Length )
{
    for(size_t i=points.size();i>0;--i)
    {
        vtx_t &r = points[i]->r;
        r.z -= 0.5;

        r.x *= Length;
        r.y *= Length;
        r.z *= Length;
    }

    for(size_t j=triangles.size();j>0;--j)
    {
        triangles[j].recompute();
    }
}


void Fish:: generateHead(  double Zmax, size_t N )
{
    clear();
    if(Zmax<=0||Zmax>0.5)
        throw exception("Invalid Zmax=%g",Zmax);

    zfunction<double> rfn( this, & Profile::getZ, Zmax);
    zfind<double>     solve(0);
    const double      rho_max =  solve(rfn.call,0,1);
    std::cerr << "rho_max=" << rho_max << std::endl;
    N = max_of<size_t>(1,N);


    const double delta = rho_max / (N+1);

    //--------------------------------------------------------------------------
    //
    // Outer Shell
    //
    //--------------------------------------------------------------------------
    // compute the slices position
    for(size_t i=1;i<=N;++i)
    {
        const double ratio = i*delta;
        const pSlice pS( new Slice( getZ(ratio) ) );
        slices.push_back(pS);
    }


    // compute the number of points per slice M = 2*n+2;
    const size_t n = max_of<size_t>(2,ceil(maxP/delta))-2;
    const size_t M = 2+2*n;
    std::cerr << "\tn=" << n << ", M=" << M << std::endl;


    // head point
    pPoint p0( new Point() );
    points.push_back(p0);

    // body points
    for(size_t i=1;i<=N;++i)
    {

        Slice &slice = *slices[i];

        const double w = W(slice.z);
        const double h = H(slice.z);

        //std::cerr << "w=" << width << ", h=" << height << std::endl;
        for(size_t j=0;j<M;++j)
        {
            const double theta = (j*numeric<double>::two_pi)/M;
            pPoint pp( new Point() );

            pp->r.x = w * cos(theta);
            pp->r.y = h * sin(theta);
            pp->r.z = slice.z;

            slice.points.push_back(pp);
            points.push_back(pp);
        }
        
    }

    std::cerr << "-- Computing Triangles" << std::endl;
    std::cerr << "\t Head..." << std::endl;

    // head
    {
        const Slice &slice = *slices[1];
        for(size_t i=1;i<=M;++i)
        {
            size_t   ip = i+1;
            if(ip>M) ip = 1;
            const Triangle tr(p0,slice.points[i],slice.points[ip]);
            triangles.push_back(tr);
        }
    }

    std::cerr << "\t Body..." << std::endl;
    // inside
    for(size_t j=1;j<N;++j)
    {
        const array<pPoint> &P0 = slices[j]->points;
        const array<pPoint> &P1 = slices[j+1]->points;

        // loop over quads
        for(size_t i=1;i<=M;++i)
        {
            size_t   ip = i+1;
            if(ip>M) ip = 1;
            const pPoint  &P00 = P0[i];
            const pPoint  &P01 = P0[ip];
            const pPoint  &P10 = P1[i];
            const pPoint  &P11 = P1[ip];

            {
                const Triangle tr(P00,P01,P11);
                triangles.push_back(tr);
            }

            {
                const Triangle tr(P00,P10,P11);
                triangles.push_back(tr);
            }
        }
    }




}



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



