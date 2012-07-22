#ifndef BUBBLE_INCLUDED
#define BUBBLE_INCLUDED 1

#include "spot.hpp"
#include "gmarker.hpp"
#include "yocto/ios/ocstream.hpp"

#if defined(HAS_MPI)
#include "yocto/mpi/mpi.hpp"
#endif

typedef ptrdiff_t BubbleID;

//! a non intersecting polygon
class Bubble :  public Point::List
{
public:
    explicit Bubble(BubbleID          who, 
                    const V2D        &Length, 
                    Real             &lambda_ref,
                    Point::Pool      &pcache, 
                    Spot::Pool       &scache,
                    GridMarker::Pool &gcache) throw();
    virtual ~Bubble() throw();
    
    const BubbleID    id;
    const PBC         pbc;         //!< from length
    const Real       &lambda;      //!< critical length
    Real              area;        //!< area, to be computed
    Real              pressure;    //!< pressure, to be broadcasted
    Spot::List        spots;    //!< keep trace of points
    GridMarker::List  markers;  //!< grid markers of inside
    bool              active;   //!< spots.size > 0 
    
    
    void   update_contour();              //!< update #points
    void   compute_values() throw();      //!< area, tangents...
    
    //! empty list and put points on circle
    void map_circle( const V2D &center, Real radius );
    
    //! empty list and put points on Cassini's shape (0<=alpha<1)
    void map_peanut( const V2D &center, Real radius, Real alpha );
    
    
    //! empty spots, mark and find out points within y_lo <= y <= y_up
    void mark_and_find_spots_within( const Real y_lo, const Real y_up );
    
    
#if defined(HAS_MPI)
    //! broadcast content from rank=0
    /**
     broadcast num_points, slaves construct bubbles
     brodacast data: pressure
     broadcast points
     */
    void dispatch( const mpi &MPI );
    
    //! assemble changed points
    void assemble( const mpi &MPI);
#endif
    
    Bubble *next;
    Bubble *prev;
    
    void save_dat( const string &filename ) const;
    void save_vtk( const string &filename ) const;
    void save_vtk_t( const string &filename ) const;
    void save_vtk_n( const string &filename ) const;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Bubble);
};




#endif
