#ifndef SPOT_INCLUDED
#define SPOT_INCLUDED 1

#include "point.hpp"
#include "yocto/core/list.hpp"

//! attach a point within a domain
struct Spot
{
    Point   *point; //!< the point
    size_t   jump;  //!< encoding: number of extra nodes to walk within the bubble
    Spot    *next;  //!< for the linked list
    Spot    *prev;  //!< for the linked list
    
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
        
        //! create a new Spot and attach p
        void append( Point *p );
        
    private:
        Pool &cache_;
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };
};



#endif