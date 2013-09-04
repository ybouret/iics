#include "yocto/fs/local-fs.hpp"
#include "yocto/math/v3d.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/string/tokenizer.hpp"
#include "yocto/exception.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/sequence/vector.hpp"

#include <iostream>

using namespace yocto;
using namespace math;

typedef double    Real;
typedef v3d<Real> V3D;

class Atom
{
public:
    V3D r;
    
    static inline bool IsSep(const char C) throw() { return C == ' ' || C == '\t'; }
    
    inline  Atom() throw() : r() {}
    inline ~Atom() throw() {}
    
    inline Atom(const Atom &atom ) throw() :
    r( atom.r )
    {
    }
    
    static inline
    Real NextReal( tokenizer &tkn, unsigned iline, const char *id )
    {
        assert(id);
        if( !tkn.get_next(IsSep) ) throw exception("%u: missing %s", iline, id);
        const string tmp = tkn.to_string();
        return strconv::to<Real>(tmp,id);
    }
    
    inline void parse( const string &line, unsigned iline)
    {
        tokenizer tkn(line);
        
        //======================================================================
        // read position
        //======================================================================
        r.x = NextReal(tkn, iline, "x");
        r.y = NextReal(tkn, iline, "y");
        r.z = NextReal(tkn, iline, "z");
        
    }
    
private:
    Atom & operator=(const Atom &);
};

typedef vector<Atom> Atoms;

class Frame : public Atoms
{
public:
    explicit Frame() throw() : Atoms()
    {
    }
    
    virtual ~Frame() throw() {}
    
    
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