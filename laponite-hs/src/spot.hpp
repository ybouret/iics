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
    
    void reset() throw();
    
    typedef cache_of<Spot> Pool;
    
    typedef cached_list<core::list_of,Spot> CoreList;
    class List : public CoreList
    {
    public:
        explicit List( Pool &cache ) throw();
        virtual ~List() throw();
       
        
        //! create a new Spot and attach p
        void attach( Point *p );
        
        
        List( const List &other );
        
    private:
        YOCTO_DISABLE_ASSIGN(List);
    };
};



#endif
