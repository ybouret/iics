#ifndef POINT_INCLUDED
#define POINT_INCLUDED 1

#include "types.hpp"
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



#endif
