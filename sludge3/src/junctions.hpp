#ifndef SLUDGE_JUNCTIONS_INCLUDED
#define SLUDGE_JUNCTIONS_INCLUDED 1

#include "grid.hpp"


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
        
        
        Junction *append(Real);
        
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
class Junctions : public Layout
{
public:
    explicit Junctions(const Grid &grid);
    virtual ~Junctions() throw();
    
    Junction::List & Vert( unit_t i ) throw();
    Junction::List & Horz( unit_t j ) throw();

private:
    size_t          jcount;
    Junction::List *jlists;
    Junction::List *jvert;   //!< width.x times
    Junction::List *jhorz;   //!< width.y times
    
};


#endif

