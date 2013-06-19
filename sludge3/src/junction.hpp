#ifndef SLUDGE_JUNCTION_INCLUDED
#define SLUDGE_JUNCTION_INCLUDED 1

#include "bubble.hpp"
#include "yocto/associative/set.hpp"

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
        const Type   type;
        const unit_t indx;
        const Real  &level;
        explicit List(Type t, const unit_t idx, const Real &v) throw();
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
    
    bool                    inside;     //!< true junction within grid
    unit_t                  lower;      //!< lower logical index on axis or INVALID
    unit_t                  upper;      //!< upper logical index on axis or INVALID
    const Bubble::Position  b_pos;      //!< default Bubble::IsInvalid, set by segmentation
    const bool              active;     //!< has an active point, depending on b_pos...
    const Tracer           *t_prev;     //!< according to bubble
    const Tracer           *t_next;     //!< according to bubble
    
    
    Junction(List &, Real, const Bubble *) throw();
    ~Junction() throw();
    
    Vertex get(void) const throw();
    Real   dist2( const Vertex &v ) const throw();
    
    
    void set_after()  const;
    void set_before() const;
    void set_active() const;
    
    class Key
    {
    public:
        const Tracer *t_prev;
        const Tracer *t_next;
        
        inline Key(const Tracer *p, const Tracer *n) throw() : t_prev(p), t_next(n) {}
        inline Key(const Junction *J) throw() : t_prev(J->t_prev), t_next(J->t_next) {}
        inline Key( const Key &k) throw() : t_prev( k.t_prev ), t_next( k.t_next ) {}
        inline ~Key() throw() {}
        inline friend
        bool operator==( const Key &lhs, const Key &rhs ) throw()
        {
            return (lhs.t_prev == rhs.t_prev) && (lhs.t_next == rhs.t_next);
        }
        
    private:
        YOCTO_DISABLE_ASSIGN(Key);
    };
    
    class Pointer
    {
    public:
        const Junction *J;
        const Key       K;
        
        inline Pointer(const Junction *j) throw() :
        J(j),
        K(J)
        {
        }
        
        inline Pointer( const Pointer &p ) throw() :
        J(p.J),
        K(p.K)
        {
        }
        
        
        inline ~Pointer() throw() {}
        
        const Key & key() const throw() { return K; }
        
    private:
        YOCTO_DISABLE_ASSIGN(Pointer);
    };
    
    typedef set<Key,Pointer> DB;
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junction);
};


#endif


