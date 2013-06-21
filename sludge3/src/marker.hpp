#ifndef SLUDGE_MARKER_INCLUDED
#define SLUDGE_MARKER_INCLUDED 1

#include "tracer.hpp"
#include "yocto/core/list.hpp"
#include "yocto/spade/array2d.hpp"

typedef array2D<Real> Array;

class Junction;

class Marker
{
public:
    static const int Tag = 4;
    Marker         *next;
    Marker         *prev;
    Tracer         *tracer;
    const size_t    shift; //!< to reconstruct tracer ID
    Real            gt;    //!< tangential pressure gradient
    Real            gn;    //!< normal gradient
    const Junction *jprev;
    const Junction *jnext;
    Vertex          v;     //!< local speed
    
    YOCTO_MAKE_OBJECT;
    Marker(Tracer *tr,const size_t s);
    ~Marker() throw();
    
    
    
    class List : public core::list_of<Marker>
    {
    public:
        explicit List() throw();
        virtual ~List() throw();
        
        void clear() throw();
        void append(Tracer *tracer, const size_t shift);
        
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Marker);
};


#endif

