#ifndef SLUDGE_JUNCTION_INCLUDED
#define SLUDGE_JUNCTION_INCLUDED 1

#include "bubble.hpp"


//! A basic junction
class Junction
{
public:
    YOCTO_MAKE_OBJECT;
    enum Type
    {
        Horz,
        Vert
    };
    
    
    //! using reference from grid !
    class List : public core::list_of<Junction>
    {
    public:
        const Type  type;
        const Real &level;
        explicit List(Type t, const Real &v) throw();
        virtual ~List() throw();
        
        const Junction * after(  const Real value ) const throw();
        const Junction * before( const Real value ) const throw();
        
        
        Junction *append(Real,const Bubble *);
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };
    
    Junction     *prev;     //!< for list
    Junction     *next;     //!< for list
    const List   &root;     //!< for level/grid
    const Real    value;    //!< location
    const Bubble *owner;    //!< owning bubble
    Real          C;        //!< average curvature
    Vertex        t;        //!< average tangent vector
    Vertex        n;        //!< average normal
    Real          pressure; //!< from curvature + owner->gamma
    
    bool              inside;     //!< true junction within grid
    unit_t            lower;      //!< lower logical index on axis or INVALID
    unit_t            upper;      //!< upper logical index on axis or INVALID
    Bubble::Position  b_pos;      //!< default Bubble::IsInvalid, set by segmentation
    
    Junction(List &, Real, const Bubble *) throw();
    ~Junction() throw();
    
    Vertex get(void) const throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junction);
};


#endif


