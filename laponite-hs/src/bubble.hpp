#ifndef BUBBLE_INCLUDED
#define BUBBLE_INCLUDED 1

#include "spot.hpp"

#if defined(HAS_MPI)
#include "yocto/mpi/mpi.hpp"
#endif

// TODO: update/update area/update differential ppties

//! a non intersecting polygon
class Bubble :  public Point::List
{
public:
    explicit Bubble( Real L, Point::Pool &pcache, Spot::Pool &scache ) throw();
    virtual ~Bubble() throw();
    
    const PBC      pbc;
    double         lambda; //!< critical length, default is 1
    double         area;   //!< area, to be broadcasted
    Spot::List     spots;  //!< keep trace of points
    bool           active; //!< spots.size > 0 
    
    
    void   update_contour();              //!< update #points
    void   compute_values() throw();      //!< area, tangents...
    
    //! empty list and put points on circle
    void map_circle( const V2D &center, Real radius );
    
    //! empty spots and find out points within y_lo <= y <= y_up
    void find_spots_within( const Real y_lo, const Real y_up );
    
    
#if defined(HAS_MPI)
    //! broadcast content from rank=0
    void dispatch( const mpi &MPI );
    
    //! assemble changed points
    void assemble( const mpi &MPI);
#endif
    
    Bubble *next;
    Bubble *prev;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
};




#endif
