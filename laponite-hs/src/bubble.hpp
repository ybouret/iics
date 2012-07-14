#ifndef BUBBLE_INCLUDED
#define BUBBLE_INCLUDED 1

#include "spot.hpp"


#if defined(HAS_MPI)
#include "yocto/mpi/mpi.hpp"
#endif

// TODO: update/update area/update differential ppties

//! a non intersecting polygon
class Bubble : public Point::List
{
public:
    explicit Bubble( Point::Pool &pcache, Spot::Pool &scache ) throw();
    virtual ~Bubble() throw();
    
    double     lambda; //!< critical length, default is 1
    double     area;   //!< area, to be broadcasted
    Spot::List spots;  //!< keep trace of points
    bool       active; //!< spots.size > 0 
    
    void   update_points(); // update #points
    void   update_area() throw(); //!< area = evaluate_area
    double evaluate_area() const throw(); //!< evaluate area, doesn't set it !
    
    //! empty list and put points on circle
    void map_circle( const V2D &center, Real radius );
    
    //! empty spots and find out points within y_lo <= y < y_up
    void build_spots( const Real y_lo, const Real y_up );
    
    
#if defined(HAS_MPI)
    //! broadcast content from rank=0
    void dispatch( mpi &MPI );
    
    //! collect changed points
    void collect(mpi &MPI);
#endif
    
    Bubble *next;
    Bubble *prev;
    
private:
    YOCTO_DISABLE_ASSIGN(Bubble);
};




#endif
