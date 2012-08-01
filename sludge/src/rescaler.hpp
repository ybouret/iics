#ifndef RESCALER_INCLUDED
#define RESCALER_INCLUDED 1

#include "bubbles.hpp"

class Rescaler
{
public:
    explicit Rescaler();
    virtual ~Rescaler() throw();
    
    

    //! upgrade on bubble
    void upgrade( Bubble &bubble );
   
    
    //! process all the bubbles !
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

    
    vector<Real> s;      //!< approximated linear abscissa
    vector<Real> ax;     //!< X coordinates
    vector<Real> ay;     //!< Y coordinates
    Real         period; //!< total s
    cache_of<abscissa>                  a_pool;
    cached_list<core::list_of,abscissa> a_list;

    //! apply pbc, compute edges, area and content to match pressure
    void build_metrics( Bubble &bubble );
    
    //! check edges length and make a list of new abscissaes
    bool need_to_refine( const Bubble &bubble );

    //! refine a precomputed metrics, if necessary
    void refine( Bubble &bubble );
       
    //! rebuild from a current metrics and a list of new abscissaes
    void rebuild( Bubble &bubble );
    
    
};



#endif
