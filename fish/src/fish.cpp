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
solv(1e-4),
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
            const double  u      = solv(zfcn.call,0,1);
            result = Compute(u);

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

Profile:: Profile( lua_State *L ) :
width(  &W, & CubiX::GetValue ),
height( &H, & CubiX::GetValue ),
W(L,"width",0),
H(L,"height",1),
zarr(NZ,0),
rmax(NZ,0)
{
    for(size_t i=1;i<=NZ;++i)
    {
        const double zz = double(i-1)/(NZ-1);
        zarr[i]  = zz;
        rmax[i]  = max_of<double>( Fabs( width(zz)), Fabs( height(zz) ) );
    }

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

Slice:: Slice(double zz,double pr) throw() :
z(zz),
perimeter(pr),
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


