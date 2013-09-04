#include "yocto/fs/local-fs.hpp"
#include "yocto/math/v3d.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/string/tokenizer.hpp"
#include "yocto/exception.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/alg/delaunay.hpp"

#include <iostream>

using namespace yocto;
using namespace math;

typedef double    Real;
typedef v3d<Real> V3D;

class Atom
{
public:
    V3D  r;
    Real voronoi;
    
    static inline bool IsSep(const char C) throw() { return C == ' ' || C == '\t'; }
    
    inline  Atom() throw() : r(), voronoi(0) {}
    inline ~Atom() throw() {}
    
    inline Atom(const Atom &atom ) throw() :
    r( atom.r ),
    voronoi( atom.voronoi )
    {
    }
    
    template <typename T>
    static inline
    T Next( tokenizer &tkn, unsigned iline, const char *id )
    {
        assert(id);
        if( !tkn.get_next(IsSep) ) throw exception("%u: missing %s", iline, id);
        const string tmp = tkn.to_string();
        return strconv::to<T>(tmp,id);
    }
    
    static inline
    void SkipNext( tokenizer &tkn, unsigned iline, const char *id )
    {
        assert(id);
        if( !tkn.get_next(IsSep) ) throw exception("%u: missing field %s", iline, id);
    }
    
    inline void parse( const string &line, unsigned iline)
    {
        tokenizer tkn(line);
        
        //======================================================================
        // read position
        //======================================================================
        r.x = Next<Real>(tkn, iline, "x");
        r.y = Next<Real>(tkn, iline, "y");
        r.z = Next<Real>(tkn, iline, "z");
        
        //======================================================================
        // read voronoi
        //======================================================================
        voronoi = Next<Real>(tkn, iline, "voronoi");
        
    }
    
private:
    Atom & operator=(const Atom &);
};

typedef vector<Atom> Atoms;
typedef vector<V3D>  Vertices;
typedef vector<iTriangle> Triangles;

class Frame : public Atoms
{
public:
    Frame *next;
    Frame *prev;
    explicit Frame() throw() :
    Atoms(),
    next(0),
    prev(0)
    {
    }
    
    void find_gas(Vertices &gas, Real vmin ) const
    {
        gas.free();
        for(size_t j=size();j>0;--j)
        {
            const Atom &atom = (*this)[j];
            if(atom.voronoi>=vmin)
                gas.push_back( atom.r );
        }
    }
    
    static
    void Triangulate( Triangles &triangles, const Vertices &gas )
    {
        delaunay<Real>::build(triangles, gas);
    }
    
    virtual ~Frame() throw() {}
    
    static Frame *LoadNext(ios::istream &fp, unsigned iline)
    {
        string line;
        if(fp.read_line(line)<0)
            return 0;
        
        
        return 0;
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Frame);
};


int main( int argc, char *argv[] )
{
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "Unhandled Exception in " << prog << std::endl;
    }
    return 1;
}