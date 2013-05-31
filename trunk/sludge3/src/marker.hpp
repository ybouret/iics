#ifndef SLUDGE_MARKER_INCLUDED
#define SLUDGE_MARKER_INCLUDED 1

#include "tracer.hpp"
#include "yocto/core/list.hpp"

class Marker
{
public:
    static const int Tag = 4;
    Marker       *next;
    Marker       *prev;
    const Tracer *tracer;
    const size_t  shift;
    
    YOCTO_MAKE_OBJECT;
    Marker(const Tracer *tr,const size_t s);
    ~Marker() throw();
    
    
    
    class List : public core::list_of<Marker>
    {
    public:
        explicit List() throw();
        virtual ~List() throw();
        
        void clear() throw();
        void append( const Tracer *tracer, const size_t shift);
        
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(List);
    };
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Marker);
};


#endif

