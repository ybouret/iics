#ifndef FISH_INCLUDED
#define FISH_INCLUDED 1

#include "yocto/math/dat/b-splines.hpp"
#include "yocto/math/v3d.hpp"
#include "yocto/counted-object.hpp"
#include "yocto/ptr/arc.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/math/fcn/zfind.hpp"



using namespace yocto;
using namespace math;

typedef CubicApproximation<double,v3d> CubiXBase;
typedef CubiXBase::vtx_t                   vtx_t;


class CubiX : public CubiXBase
{
public:
    explicit CubiX(lua_State    *L,
                   const string &table_name,
                   const size_t  iCoord);
    virtual ~CubiX() throw();

    double  GetValue(const double z); //! X or Y

private:
    zfind<double>     solv; //!< solver
    zfunction<double> zfcn; //!< for inversion
    const size_t      indx; //!< 0 -> X, 1->Y
    YOCTO_DISABLE_COPY_AND_ASSIGN(CubiX);

};

class Profile : public object
{
public:
    numeric<double>::function width;
    numeric<double>::function height;

    explicit Profile(lua_State *L );
    virtual ~Profile() throw();

    
    static const size_t NZ = 1001;

private:
    CubiX W;
    CubiX H;
public:
    vector<double> zarr; //!< support points
    vector<double> rmax; //!< max radius at this point
    YOCTO_DISABLE_COPY_AND_ASSIGN(Profile);
};




class Point : public counted_object
{
public:
    const size_t i;
    vtx_t        r;
    explicit Point() throw();
    virtual ~Point() throw();

private:
    static size_t Index;
};

typedef arc_ptr<Point> pPoint;


class Slice : public counted_object
{
public:
    const double   z;
    const double   perimeter;
    vector<pPoint> points;

    explicit Slice(double zz,double pr) throw();
    virtual ~Slice() throw();

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Slice);
};

typedef arc_ptr<Slice> pSlice;

// Triangle with normal out w.r.t [0 0 0.5]
class Triangle
{
public:
    pPoint a,b,c;
    vtx_t  G;
    vtx_t  n;
    double S;

    Triangle(const pPoint &A,
             const pPoint &B,
             const pPoint &C) throw();
    
    Triangle(const Triangle &other) throw();
    ~Triangle() throw();

    //! after a scaling...
    void recompute() throw();
    
private:
    YOCTO_DISABLE_ASSIGN(Triangle);
};



#endif

