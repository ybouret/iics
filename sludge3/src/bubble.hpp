#ifndef SLUDGE_BUBBLE_INCLUDED
#define SLUDGE_BUBBLE_INCLUDED 1

#include "marker.hpp"
#include "yocto/ios/ostream.hpp"
#include "yocto/string.hpp"


class Bubble : public Tracer::Ring
{
public:
    static const int Tag;
    
    enum Position
    {
        IsInvalid,
        IsBefore,
        IsAfter
    };
    
    
    Bubble *prev;
    Bubble *next;
    explicit Bubble( Real &lam, Real &gam, size_t uid) throw();
    virtual ~Bubble() throw();
    
    const Real         &lambda;   //!< maximum length between two vertices
    const Real         &gamma;    //!< surface tension
    Vertex              G;        //!< barycenter    : +2 Real
    Real                area;     //!< area          : +1 Real
    Real                pressure; //!< pressure      : +1 Real, default is 1
    static const size_t NumReals = 4;
    
    size_t              flags;
    const size_t        UID;      //!< for segmentation
    Marker::List        markers;  //!< for MPI processing
    
    void clear() throw();
    
    void collect_markers( const Real ymin, const Real ymax);
    
    
    void save_dat( const string &fn ) const;
    void save_t( const string &fn ) const; //!< tangents
    void save_n( const string &fn ) const; //!< curv + normal
    void save_p( const string &fn ) const; //!< probing
    void save_all( const string &pfx ) const;
    
    
    //! more than three points !
    /**
     compute edges and distances, update area.
     */
    void init_contour() throw();
    
    
    //! must have an init_contour before
    void adjust_contour();
    
    
    //! compute tangent/normal/curvature after an [init|auto]_contour
    void compute_curvatures();

    
    //! auto_contour() + compute_curvatures()
    void regularize();
    
    //! used for MPI debugging
    void hash_bubble( Hasher &h ) const throw();
    
    Tracer *append(); //!< helper
    
    void    append( const Vertex v ); //!< use a copy or vertex
    
    
    void dispatch_junctions();
    

private:

    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
    
};

#endif

