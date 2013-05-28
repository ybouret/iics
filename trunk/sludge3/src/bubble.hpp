#ifndef SLUDGE_BUBBLE_INCLUDED
#define SLUDGE_BUBBLE_INCLUDED 1

#include "tracer.hpp"
#include "yocto/ios/ostream.hpp"
#include "yocto/string.hpp"

class Bubble : public Tracer::Ring
{
public:
    static const int Tag;
    
    Bubble *prev;
    Bubble *next;
    explicit Bubble( Real &lam ) throw();
    virtual ~Bubble() throw();
    
    const Real &lambda;   //!< maximum length between two vertices
    Vertex      G;        //!< barycenter    : +2 Real
    Real        area;     //!< area          : +1 Real
    Real        pressure; //!< pressure      : +1 Real, default is 1
    static const size_t NumReals = 4;
    
    void save_dat( ios::ostream &fp ) const;
    void save_dat( const string &fn ) const;
    
    //! more than three points !
    /**
     compute edges and distances, update area.
     */
    void init_contour() throw();
    
    //! once init_contour is ok, respect lambda and init new contour
    void auto_contour();
    
    void hash_bubble( Hasher &h ) const throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
    
};

#endif

