#ifndef BUBBLE_INCLUDED
#define BUBBLE_INCLUDED 1

#include "types.hpp"
#include "yocto/core/list.hpp"
#include "yocto/core/pool.hpp"
#include "yocto/core/circular-list.hpp"

class Point : public V2D
{
public:
    Real     s_next;  //!< distance to next point
    
    explicit Point() throw(); //!< x=y=0, domain=-1, d2next=0
    virtual ~Point() throw(); //!< do nothing
    
    Point( const Point & ) throw(); //! copy
    Point & operator=( const Point & ) throw();
    
    Point *next; //!< for linked list
    Point *prev; //!< for linked list
    
    typedef core::pool_of<Point> CorePool;
    class Pool : public CorePool
    {
    public:
        explicit Pool() throw();
        explicit Pool(size_t cache_size);
        virtual ~Pool() throw();
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Pool);
    };
    
    typedef core::clist_of<Point> CoreList;
    class List : public CoreList
    {
    public:
        explicit List( Pool &cache ) throw(); //!< empty list with a new cache
        virtual ~List() throw();              //!< put back into cache
        
        void empty() throw(); //!< empty the points back into cache
        
        Point *create();  //!< get a new point (from cache if possible)
        void   remove( Point *p ) throw(); //!< unlink and cache the point
        
    private:
        Pool &cache_;
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };
        
};

struct Spot
{
    Point *point;
    Spot  *next;
    Spot  *prev;
    
    typedef core::pool_of<Spot> CorePool;
    class Pool : public CorePool
    {
    public:
        explicit Pool() throw();
        virtual ~Pool() throw();
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Pool);
    };
    
    typedef core::list_of<Spot> CoreList;
    class List : public CoreList
    {
    public:
        explicit List( Pool &cache ) throw();
        virtual ~List() throw();
        void     empty() throw();
        
        void append( Point *p );
        
    private:
        Pool &cache_;
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };
};


//! a non intersecting polygon
class Bubble : public Point::List
{
public:
    explicit Bubble( Point::Pool &cache, Spot::Pool &spot_cache ) throw();
    virtual ~Bubble() throw();
    
    double     lambda; //!< critical length, default is 1
    double     area;   //!< area
    Spot::List spots;  //!< keep trace of points
    
    void   update();
    double evaluate_area() const throw(); //!< evaluate area, doesn't set it !
    
    //! empty list and put points on circle
    void map_circle( const V2D &center, Real radius );
    
    //! empty spots and find out point within y_lo <= y < y_up
    void build_spots( const Real y_lo, const Real y_up );
    
private:
    YOCTO_DISABLE_ASSIGN(Bubble);
};




#endif
