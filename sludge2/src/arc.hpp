#ifndef ARC_INCLUDED
#define ARC_INCLUDED 1

#include "./types.hpp"
#include "yocto/math/kernel/linsys.hpp"


class Arc
{
public:
    Arc() throw();
    ~Arc() throw();
    
    mutable Vertex a;
    mutable Vertex b;
    mutable Vertex c;
    
    Vertex operator()( const Real mu ) const throw();
    
    
    Vertex r0;
    Vertex t0;
    Real   C0;
    
    Vertex r1;
    Vertex t1;
    Real   C1;
    
    mutable Vertex dr;
    
    // extract a,b,c from U
    void load( const array<Real> &U ) const;
    
    //! fill rows 3 and 4, and compute dr, and initialize U
    void fill( matrix<Real> &P, array<Real> &U) const;
    
    //! fill rows 5 and 6 of P,
    void estimate( matrix<Real> &P, array<Real> &F, const array<Real> &U ) const;
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Arc);
};

class ArcSolver
{
public:
    ArcSolver();
    ~ArcSolver() throw();
    
    void compute( const Arc &arc );
    
private:
    matrix<Real> P;
    matrix<Real> iP;
    vector<Real> U;
    vector<Real> F;
    vector<Real> h;
    linsys<Real> ls;
    
    YOCTO_DISABLE_COPY_AND_ASSIGN(ArcSolver);
};

#endif

