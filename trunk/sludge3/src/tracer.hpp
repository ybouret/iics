#ifndef SLUDGE_TRACER_INCLUDED
#define SLUDGE_TRACER_INCLUDED 1

#include "types.hpp"
#include "yocto/core/clist.hpp"

class Tracer : public object
{
public:
    static const int Tag;
    
    Tracer *prev;
    Tracer *next;
    
    Vertex pos;  //!< current pos +2 Real
    Vertex edge; //!< to next     +2 Real
    Real   dist; //!< to next     +1 Real
    
    static const size_t NumReals = 5;
    
    explicit Tracer() throw();
    explicit Tracer( const Vertex v ) throw();
    virtual ~Tracer() throw();
    
    void hash_tracer( Hasher &h ) const throw();
    
    class Ring : public core::clist_of<Tracer>
    {
    public:
        explicit Ring() throw();
        virtual ~Ring() throw();
        void hash_ring( Hasher &h) const throw();
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Ring);
    };
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Tracer);
};



#endif
