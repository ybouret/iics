#ifndef ARC_INCLUDED
#define ARC_INCLUDED 1

#include "./junction.hpp"

#if 0
class ArcPoint
{
public:
    const Vertex r;     //!< position
    const Vertex t;     //!< tangent vector
    const Real   P;     //!< pressure
    const Real   alpha; //!< tangential gradient
    const Real   beta;  //!< orthogonal gradient
    const Real   theta; //!< angle ot t
    mutable Real delta; //!< angle difference, computed later
    mutable Real C0;
    mutable Real S0;
    mutable Real C1;
    mutable Real S1;
    mutable Real C2;
    mutable Real S2;
    
    ArcPoint(const Vertex &a_r,
             const Vertex &a_t,
             const Real    a_P,
             const Real    gt,
             const Real    gn);
    
    ~ArcPoint() throw();
    
    void computeIntegrals() const;
    
private:
    Real dC0( Real ) const; // cos(theta+mu*delta)
    Real dS0( Real ) const; // sin(theta+mu*delta)
    Real dC1( Real ) const; // mu*cos(theta+mu*delta)
    Real dS1( Real ) const; // mu*sin(theta+mu*delta)
    Real dC2( Real ) const; // mu^2*cos(theta+mu*delta)
    Real dS2( Real ) const; // mu^2*sin(theta+mu*delta)
    YOCTO_DISABLE_COPY_AND_ASSIGN(ArcPoint);
};

#include "yocto/math/kernel/svd.hpp"


class Arc
{
public:
    Arc(const ArcPoint &a_A,
        const ArcPoint &a_Q,
        const ArcPoint &a_B) throw();
    
    
    ~Arc() throw();
    
    const ArcPoint &A;
    const ArcPoint &Q;
    const ArcPoint &B;
    
    const Vertex AQ;
    const Vertex BQ;
    
    Real I0;
    Real J0;
    Real I1;
    Real J1;
    Real I2;
    Real J2;
    
    Real I0p;
    Real J0p;
    Real I1p;
    Real J1p;
    Real I2p;
    Real J2p;
    
    Real eta;
    Real etap;
    
    void load( matrix<Real> &H, array<Real> &U, matrix<Real> &JK ) const throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Arc);
};


class ArcSolver
{
public:
    static const size_t N=5; //! lagrange multipliers
    static const size_t M=8; //!< unknown parameters
    ArcSolver();
    ~ArcSolver() throw();
    
    void operator()( const Arc &arc );
    
private:
    matrix<Real> H;
    vector<Real> U;
    vector<Real> L;
    vector<Real> W;
    matrix<Real> V;
    matrix<Real> JK;
    vector<Real> X;
    YOCTO_DISABLE_COPY_AND_ASSIGN(ArcSolver);
    
};
#endif


#endif

