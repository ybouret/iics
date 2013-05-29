#ifndef SLUDGE_JUNCTIONS_INCLUDED
#define SLUDGE_JUNCTIONS_INCLUDED 1

#include "grid.hpp"
#include "bubbles.hpp"


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
        
        
        void append(Real);
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };

    Junction   *prev;
    Junction   *next;
    const List &root;
    const Real  value;
    
    explicit Junction(List &, Real) throw();
    virtual ~Junction() throw();
    
        
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junction);
};


//! handle Vert/Horz junctions from grid
class Junctions 
{
public:
    explicit Junctions(Grid &g);
    virtual ~Junctions() throw();
    
    Junction::List & Vert( unit_t i ) throw();
    Junction::List & Horz( unit_t j ) throw();

    size_t load( const Bubble &bubble ); //!< append new junctions, return status
    void   clear() throw();              //!< clear all junctions
    void   sort();
    
    
    
    void save_dat( const string &fn) const;
    
    const Grid     &grid;
    const size_t    num_lists;
   

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junctions);
    size_t          jcount;
    Junction::List *jlists;
    Junction::List *jvert;   //!< width.x times
    Junction::List *jhorz;   //!< width.y times
    
    void __load( const Bubble &bubble, const Tracer *u);
    void __loadJ( const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q, const Coord &Q);
    
    void __loadHorz(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q);
    void __loadVert(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q);

};


#endif

