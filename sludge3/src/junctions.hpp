#ifndef SLUDGE_JUNCTIONS_INCLUDED
#define SLUDGE_JUNCTIONS_INCLUDED 1

#include "grid.hpp"
#include "bubbles.hpp"
#include "yocto/spade/array2d.hpp"

typedef array2D<Real> Array;

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
    
    bool   inside;
    unit_t lower;
    unit_t upper;
    
    Junction(List &, Real, const Bubble *) throw();
    ~Junction() throw();
    
    Vertex get(void) const throw();
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junction);
};


//! handle Vert/Horz junctions from grid
class Junctions
{
public:
    const Grid     &grid;
    const size_t    num_lists;
    
    explicit Junctions(Grid &g);
    virtual ~Junctions() throw();
    
    Junction::List & Vert( unit_t i ) throw();
    Junction::List & Horz( unit_t j ) throw();
    
    const Junction::List & Vert( unit_t i ) const throw();
    const Junction::List & Horz( unit_t j ) const throw();
    
    void   inter(Bubble &bubble);  //!< create geometrical intersections, set bubble.flags
    void   clear() throw();        //!< clear all junctions
    void   sort();                 //!< sort for segmentation
    
    
    void load( Bubbles &bubbles );
    
    
    void save_dat( const string &fn) const;
    void save_t( const string &fn ) const;
    void save_n( const string &fn ) const;
    
    
    void segment( Array &B ) const; //!< fill array with owner
    
    void save_inside( const Array &B, const string &fn ) const;
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junctions);
    size_t          jcount;
    Junction::List *jlists;
    Junction::List *jvert;   //!< width.x times
    Junction::List *jhorz;   //!< width.y times
    
    
    void __intersect(const Bubble &bubble, const Tracer *u);
    
    // create junction with average pressure
    Junction *__interHorz(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q, Real &alpha);
    Junction *__interVert(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q, Real &alpha);
    
    void __updateJunction( Junction *J, const Real alpha, const Tracer *u, const Tracer *v);
    
    
};


#endif

