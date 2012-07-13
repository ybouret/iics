#ifndef BUBBLE_INCLUDED
#define BUBBLE_INCLUDED 1

#include "types.hpp"
#include "yocto/core/list.hpp"
#include "yocto/core/pool.hpp"

class Point : public object, public V2D
{
public:
    int      domain;
    explicit Point() throw(); //!< x=y=0, domain=-1
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
    
    typedef core::list_of<Point> CoreList;
    class List : public CoreList
    {
    public:
        explicit List( Pool &cache ) throw(); //!< empty list with a new cache
        virtual ~List() throw();              //!< put back into cache
        
        Point *New();  //!< get a new point (from cache if possible)
        
    private:
        Pool &cache_;
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };
    
};


class Polygon : public Point::List
{
public:
    explicit Polygon( Point::Pool &cache ) throw();
    virtual ~Polygon() throw();
    
    double area;
    void   Update() throw();
    
    
private:
    YOCTO_DISABLE_ASSIGN(Polygon);
};




#endif
