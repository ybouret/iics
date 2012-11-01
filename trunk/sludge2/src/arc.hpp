#ifndef ARC_INCLUDED
#define ARC_INCLUDED 1

#include "./types.hpp"
#include "yocto/math/kernel/matrix.hpp"

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
    
    mutable Vertex delta_r;
    mutable Real   delta;
    mutable Real   delta2;
    
    //! extract a,b,c,alpha,beta from U
    void load( const array<Real> &U ) const;
    
    //! compute delta_r and init U
    void init( array<Real> &U ) const;
    
    
    void func( array<Real>  &F, const array<Real> &U ) const;
    void fjac( matrix<Real> &J, const array<Real> &U ) const;
    
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
    vector<Real> U;
    
    YOCTO_DISABLE_COPY_AND_ASSIGN(ArcSolver);
};

#endif

