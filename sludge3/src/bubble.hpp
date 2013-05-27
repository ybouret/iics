#ifndef SLUDGE_BUBBLE_INCLUDED
#define SLUDGE_BUBBLE_INCLUDED 1

#include "tracer.hpp"
#include "yocto/ios/ostream.hpp"
#include "yocto/string.hpp"

class Bubble : public Tracer::Ring
{
public:
    explicit Bubble( Real &lam ) throw();
    virtual ~Bubble() throw();
    
    const Real &lambda; //!< maximum length between two vertices
    
    void save_dat( ios::ostream &fp ) const;
    void save_dat( const string &fn ) const;
    
    //! more than three points !
    void init_contour() throw();
    
    //! once init_contour is ok, respect lambda
    void auto_contour();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
    
};

#endif

