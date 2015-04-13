#ifndef FISH_INCLUDED
#define FISH_INCLUDED 1

#include "yocto/math/dat/b-splines.hpp"
#include "yocto/math/v3d.hpp"
#include "yocto/counted-object.hpp"
#include "yocto/ptr/arc.hpp"

using namespace yocto;
using namespace math;

typedef CubicApproximation<double,v3d> CubiX;
typedef CubiX::vtx_t                   vtx_t;

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

