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


//==============================================================================
//
//! Cubic BSplines Approximation
//
//==============================================================================
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

//==============================================================================
//
//! Profile Metrics
//
//==============================================================================
class Profile : public object
{
public:
    numeric<double>::function W;
    numeric<double>::function H;

    explicit Profile(lua_State *L );
    virtual ~Profile() throw();

    
    static const size_t NZ = 1001;

    //! required Z position for a given ratio of arc
    double getZ( const double ratio );

    //! compute the perimeter
    double computePerimeter( const double z );



private:
    CubiX width;
    CubiX height;
    double minusPerimeter( const double );
    
public:
    vector<double> zarr;   //!< support points
    vector<double> rmax;   //!< max radius at this point
    vector<double> arcL;   //!< arc length, 0->1
    const double   maxL;   //!< max arc length
    const double   maxP;   //!< max perimeter
    YOCTO_DISABLE_COPY_AND_ASSIGN(Profile);
};


//==============================================================================
//
//! a data point
//
//==============================================================================
class Point : public counted_object
{
public:
    const size_t i;
    vtx_t        r;
    explicit Point() throw();
    virtual ~Point() throw();

    static size_t Index;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Point);
};

typedef arc_ptr<Point> pPoint;


//==============================================================================
//
//! an elliptical slice
//
//==============================================================================
class Slice : public counted_object
{
public:
    const double   z;
    vector<pPoint> points;

    explicit Slice(double zz) throw();
    virtual ~Slice() throw();

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Slice);
};

typedef arc_ptr<Slice> pSlice;


//==============================================================================
//
//! Triangle with normal out w.r.t [0 0 0.5]
//
//==============================================================================
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

    //! change orientation
    void inverse() throw();

private:
    YOCTO_DISABLE_ASSIGN(Triangle);
};


//==============================================================================
//
//!
//
//==============================================================================
class Fish : public Profile
{
public:
    explicit Fish( lua_State *L );
    virtual ~Fish() throw();

    vector<pSlice>   slices;
    vector<pPoint>   points;
    vector<Triangle> triangles;

    void clear() throw();
    


    //! generate shell between a and b
    void generateShell( size_t N );

    //! rescale all, assuming initial z\in[0:1]
    void centerAndRescaleBy( double Length );

    void generateHead( double Zmax, size_t N, double thickness );

    void generateTail( double Zmax, size_t N );


    void save_vtk( const string &filename ) const;
    void save_stl( const string &filename ) const;

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Fish);
};


#endif

