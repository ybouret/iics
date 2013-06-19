#ifndef SLUDGE_JUNCTIONS_INCLUDED
#define SLUDGE_JUNCTIONS_INCLUDED 1

#include "grid.hpp"
#include "bubbles.hpp"
#include "junction.hpp"
#include "yocto/spade/array2d.hpp"



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
    void save_all( const string &pfx) const;
    
    
    
    void segment( Array &B ) const; //!< fill array with owner
    
    void save_inside_of( const Array &B, const string &fn ) const;
    
    size_t count_all() const throw();
    void   to_curve( array<Real> &cx, array<Real> &cy ) const throw();
    
    
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Junctions);
    size_t          jcount;
    Junction::List *jlists;
    Junction::List *jvert;   //!< width.x times
    Junction::List *jhorz;   //!< width.y times
    
    void __intersect(Bubble &bubble, const Tracer *u);
    
    // create junction with average pressure/curvature
    Junction *__interHorz(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q, Real &alpha);
    
    // create junction with average pressure/curvature
    Junction *__interVert(const Bubble &bubble, const Vertex &p, const Coord &P, const Vertex &q, Real &alpha);
    
    void __updateJunction( Junction *J, const Real alpha, const Tracer *u, const Tracer *v);
    
    
};


#endif

