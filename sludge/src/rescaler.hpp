#ifndef RESCALER_INCLUDED
#define RESCALER_INCLUDED 1

#include "bubbles.hpp"
#include "yocto/math/dat/trigonometric.hpp"


class Rescaler
{
public:
    explicit Rescaler();
    virtual ~Rescaler() throw();
    
    
    void process( Bubbles &bubbles );
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Rescaler);
    
    struct abscissa
    {
        abscissa *next;
        abscissa *prev;
        Real      s;
        void reset() throw();
    };

    
    linsys<Real> solver;
    vector<Real> s;      //!< approximated linear abscissa
    vector<Real> ax;     //!< X coordinates
    vector<Real> ay;     //!< Y coordinates
    Real         period; //!< total s
    vector<Real> theta;  //!< corresponding periodic abscissa

    //! refine a precomputed metrics, if necessary
    void refine( Bubble &bubble );
    
    cache_of<abscissa>                  a_pool;
    cached_list<core::list_of,abscissa> a_list;
    
    void build_metrics( Bubble &bubble );
    void rebuild( Bubble &bubble );
    
};



#endif
