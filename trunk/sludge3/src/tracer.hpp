#ifndef SLUDGE_TRACER_INCLUDED
#define SLUDGE_TRACER_INCLUDED 1

#include "types.hpp"
#include "yocto/spade/in2d.hpp"
#include "yocto/core/clist.hpp"

using namespace spade;

#define SLUDGE_INVALID_COORD (limit_of<unit_t>::maximum)


class Tracer
{
public:
    static const int Tag;
    
    YOCTO_MAKE_OBJECT;
    
    Tracer *prev;
    Tracer *next;
    
    Vertex pos;  //!< current pos   : +2 Real
    Vertex edge; //!< to next       : +2 Real
    Real   dist; //!< to next       : +1 Real
    Vertex t;    //!< tangent vector: +2 Real
    Vertex n;    //!< normal  vector: +2 real
    Real   C;    //!< curvature     : +1 real
    static const size_t NumReals = 10;
    
    explicit Tracer() throw();
    explicit Tracer( const Vertex &v ) throw();
    ~Tracer() throw();
    
    void hash_tracer( Hasher &h ) const throw();
    
    void compute_order1(); //!< tau and n, keep |dM| as speed
    void compute_order2(); //!< evaluate curvature
    
    
    
    class Ring : public core::clist_of<Tracer>
    {
    public:
        explicit Ring() throw();
        virtual ~Ring() throw();
        void hash_ring( Hasher &h) const throw();
        Real __area() const throw();
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Ring);
    };
    Real            speed;   //!< local dM/dt store
    mutable coord2D coord;   //!< on a grid
    mutable size_t  flags;   //!< position flag
    
private:
    
    YOCTO_DISABLE_COPY_AND_ASSIGN(Tracer);
};



#endif
