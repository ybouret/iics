#ifndef SLUDGE_JUNCTIONS_INCLUDED
#define SLUDGE_JUNCTIONS_INCLUDED 1

#include "grid.hpp"
#include "bubbles.hpp"
#include "yocto/spade/array2d.hpp"

typedef array2D<Real> Array;

//! A basic junction
class Junction : public object
{
public:
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
        
        
        void append(Real,const Bubble *);
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };

    Junction     *prev;
    Junction     *next;
    const List   &root;
    const Real    value;
    const Bubble *owner;
    
    explicit Junction(List &, Real, const Bubble *) throw();
    virtual ~Junction() throw();
    
    Coord operator()(void) const throw();
    
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

    void   inter(Bubble &bubble);  //!< create geometrical intersections, set bubble.flags
    void   clear() throw();        //!< clear all junctions
    void   sort();                 //!< sort for segmentation
    
    
    void load( Bubbles &bubbles );
    
    
    void save_dat( const string &fn) const;
    
    
    
    void segment( Array &B ) const; //!< fill array with owner

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junctions);
    size_t          jcount;
    Junction::List *jlists;
    Junction::List *jvert;   //!< width.x times
    Junction::List *jhorz;   //!< width.y times
    
    
    void __intersect(const Bubble &bubble, const Tracer *u);
    void __interHorz(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q);
    void __interVert(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q);

};


#endif

