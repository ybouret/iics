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
    explicit Bubble( Real &lam, size_t uid) throw();
    virtual ~Bubble() throw();
    
    const Real         &lambda;   //!< maximum length between two vertices
    Vertex              G;        //!< barycenter    : +2 Real
    Real                area;     //!< area          : +1 Real
    Real                pressure; //!< pressure      : +1 Real, default is 1
    static const size_t NumReals = 4;
    size_t              flags;
    const size_t        UID;      //!< for segmentation
    
    void save_dat( const string &fn ) const;
    void save_t( const string &fn ) const;
    void save_n( const string &fn ) const;
    
    //! more than three points !
    /**
     compute edges and distances, update area.
     */
    void init_contour() throw();
    
    //! once init_contour is ok, respect lambda and init new contour
    void auto_contour();
    
    //! compute tangent/normal/curvature after an [init|auto]_contour
    void compute_curvatures();

    
    //! used for MPI debugging
    void hash_bubble( Hasher &h ) const throw();
    
    Tracer *append();
    void    append( const Vertex v ); //!< use a copy or vertex
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
    
};

#endif

