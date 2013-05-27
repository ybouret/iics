#ifndef SLUDGE_TRACER_INCLUDED
#define SLUDGE_TRACER_INCLUDED 1

#include "types.hpp"
#include "yocto/core/clist.hpp"

class Tracer : public object
{
public:
    Tracer *prev;
    Tracer *next;
    
    Vertex pos;  //!< current pos
    Vertex edge; //!< to next
    Real   dist; //!< to next
    
    
    explicit Tracer() throw();
    explicit Tracer( const Vertex v ) throw();
    virtual ~Tracer() throw();
    
    void hash( Hasher &h ) const throw();
    
    class Ring : public core::clist_of<Tracer>
    {
    public:
        explicit Ring() throw();
        virtual ~Ring() throw();
        void hash( Hasher &h) const throw();
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Ring);
    };
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Tracer);
};



#endif
