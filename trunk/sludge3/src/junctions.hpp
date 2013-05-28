#ifndef SLUDGE_JUNCTIONS_INCLUDED
#define SLUDGE_JUNCTIONS_INCLUDED 1

#include "grid.hpp"


class Junction : public object
{
public:
    enum Type
    {
        Horz,
        Vert
    };
    
    class List : public core::list_of<Junction>
    {
    public:
        const Type  type;
        const Real &level;
        explicit List(Type t, const Real &v) throw();
        virtual ~List() throw();
        
        
        Junction *append();
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };

    Junction   *prev;
    Junction   *next;
    const List &root;
    
    explicit Junction(List &) throw();
    virtual ~Junction() throw();
    
        
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junction);
};

class Junctions : public Layout
{
public:
    explicit Junctions(const Grid &grid);
    virtual ~Junctions() throw();
    
private:
    size_t          jcount;
    Junction::List *jlists;
    Junction::List *jvert;   //!< width.x times
    Junction::List *jhorz;   //!< width.y times
    
};


#endif

