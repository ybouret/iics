#ifndef ARC_INCLUDED
#define ARC_INCLUDED 1

#if 0
#include "./junction.hpp"

class ArcPoint
{
public:
    const Vertex r;     //!< position
    const Vertex t;     //!< tangent vector
    const Real   P;     //!< pressure
    const Real   alpha; //!< tangential gradient
    const Real   beta;  //!< orthogonal gradient
    const Real   theta; //!< angle of t
    mutable Real delta; //!< angle difference, computed later
    
    
    ArcPoint(const Vertex &a_r,
             const Vertex &a_t,
             const Real    a_P,
             const Real    gt,
             const Real    gn);
    
    ~ArcPoint() throw();
    
    void computeIntegrals() const;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(ArcPoint);
};

#include "yocto/math/kernel/svd.hpp"
#include "yocto/math/kernel/linsys.hpp"


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
    const Vertex QB;
    const Vertex AB;
    const Real   lambda;
    const Real   omlam;  //!< 1-lambda
    const Real   lambda2;
    const Real   lambda3;
    Real theta(Real mu) const throw();
    Real   theta_AQ(Real mu) const throw();
    Real   theta_QB(Real mu) const throw();
    Vertex rdot_AQ(Real mu) const throw();
    Vertex rdot_QB(Real mu) const throw();
    
    void load( matrix<Real> &K, array<Real> &Z) const throw();
    mutable size_t nu;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Arc);
    Real dI_AQ(Real) const throw();
    Real dJ_AQ(Real) const throw();
    
    Real dI_QB(Real) const throw();
    Real dJ_QB(Real) const throw();
};


class ArcSolver
{
public:
    static const size_t NC=5;
    static const size_t NX=6;
    explicit ArcSolver();
    virtual ~ArcSolver() throw();
    
    void operator()( const Arc &arc );
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(ArcSolver);
    
    matrix<Real> K;
    matrix<Real> U; //!< K' => SVD_U
    vector<Real> W; //!< SVD_W
    matrix<Real> V; //!< SVD_V
    vector<Real> Z;
    matrix<Real> J;
    matrix<Real> JU;
    matrix<Real> H;
    //vector<Real> iW; //!< approx 1/W
    vector<Real> Y;  //!< for computation
    linsys<Real> solve;
    vector<Real> X;
    
    Real dI_AQ(Real) const throw();
    Real dJ_AQ(Real) const throw();
    
};

#endif
#endif



