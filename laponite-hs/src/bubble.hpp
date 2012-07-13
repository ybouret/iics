#ifndef BUBBLE_INCLUDED
#define BUBBLE_INCLUDED 1

#include "types.hpp"
#include "yocto/core/list.hpp"
#include "yocto/core/pool.hpp"
#include "yocto/core/circular-list.hpp"

class Point : public V2D
{
public:
    int      domain; //!< MPI style owner
    Real     d2next; //!< distance to next point
    
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
        
    private:
        Pool &cache_;
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };
    
};

typedef Point::Pool PCache;

//! a non intersectiong polygon
class Bubble : public Point::List
{
public:
    explicit Bubble( Point::Pool &cache ) throw();
    virtual ~Bubble() throw();
    
    double lambda; //!< critical length, default is 1
    double area;   //!< area
    
    void   update();
    
    //! empty list and put points on circle
    void map_circle( const V2D &center, Real radius );
    
private:
    YOCTO_DISABLE_ASSIGN(Bubble);
};




#endif
