#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/v3d.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/string/tokenizer.hpp"
#include "yocto/exception.hpp"
#include "yocto/math/alg/delaunay.hpp"
#include "yocto/ordered/sorted-vector.hpp"
#include "yocto/sort/heap.hpp"
#include "yocto/sys/wtime.hpp"

using namespace yocto;
using namespace math;

////////////////////////////////////////////////////////////////////////////////
//
// Types definition
//
////////////////////////////////////////////////////////////////////////////////
typedef double                 Real;
typedef v3d<Real>              V3D;
typedef vector<V3D>            Vertices;
typedef vector<iTriangle>      Triangles;
typedef sorted_vector<size_t>  Sorted;



////////////////////////////////////////////////////////////////////////////////
//
// One Atom
//
////////////////////////////////////////////////////////////////////////////////
class Atom
{
public:
    unsigned id;
    V3D      r;
    V3D      v;
    Real     voronoi;
    
    static inline bool IsSep(const char C) throw() { return C == ' ' || C == '\t'; }
    
    inline  Atom() throw() : id(0), r(), v(), voronoi(0) {}
    inline ~Atom() throw() {}
    
    inline Atom(const Atom &atom ) throw() :
    id( atom.id ),
    r( atom.r ),
    v( atom.v ),
    voronoi( atom.voronoi )
    {
    }
    
    static inline
    int Compare( const Atom &lhs, const Atom &rhs ) throw()
    {
        return ptrdiff_t(lhs.id) - ptrdiff_t(rhs.id);
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
    
    inline void parse_full( const string &line, unsigned iline)
    {
        tokenizer tkn(line);
        
        //======================================================================
        // read id
        //======================================================================
        id = Next<size_t>(tkn,iline,"id");
        
        //======================================================================
        // read position
        //======================================================================
        r.x = Next<Real>(tkn, iline, "x");
        r.y = Next<Real>(tkn, iline, "y");
        r.z = Next<Real>(tkn, iline, "z");
        
        //======================================================================
        // read velocity
        //======================================================================
        v.x = Next<Real>(tkn, iline, "vx");
        v.y = Next<Real>(tkn, iline, "vy");
        v.z = Next<Real>(tkn, iline, "vz");
        
        
        //======================================================================
        // read voronoi
        //======================================================================
        voronoi = Next<Real>(tkn, iline, "voronoi");
    }
    
private:
    Atom & operator=(const Atom &);
};

YOCTO_SUPPORT_NO_DESTRUCT(Atom)

typedef vector<Atom>           Atoms;

////////////////////////////////////////////////////////////////////////////////
//
// One Frame
//
////////////////////////////////////////////////////////////////////////////////
enum LoadMode
{
    LoadLiquid,
    LoadMembrane
};

class Frame : public Atoms
{
public:
    Frame *next;
    Frame *prev;
    Real   runtime;
    Real   bps;
    wtime  chrono;
    
    explicit Frame() throw() :
    Atoms(),
    next(0),
    prev(0),
    runtime(0),
    bps(0),
    chrono()
    {
        chrono.start();
    }
    
    // Handle Memory
    void prepare(size_t na)
    {
        free();
        ensure(na);
        bps   = 0;
    }
    
    void find_gas(Vertices &gas, Real vmin, Sorted &gid ) const
    {
        gas.free();
        gid.free();
        for(size_t j=size();j>0;--j)
        {
            const Atom &atom = (*this)[j];
            if(atom.voronoi>=vmin)
            {
                gas.push_back( atom.r );
                gas.back().z = 0;
                if(!gid.insert(j))
                    throw exception("internal indices failure");
            }
        }
    }
    
    static
    void Triangulate( Triangles &triangles, const Vertices &gas )
    {
        delaunay<Real>::build(triangles, gas);
    }
    
    virtual ~Frame() throw() {}
    
    static bool ReadLine( string &line, ios::istream &fp, unsigned &iline )
    {
        line.clear();
        if( fp.read_line(line) < 0 )
            return false;
        
        ++iline;
        return true;
    }
    
    
    bool load_next(ios::icstream &fp,
                   unsigned      &iline,
                   LoadMode       mode)
    {
        string          line;
        const Real      t_ini = chrono.query();
        const ptrdiff_t p_ini = fp.tell();
        
        // read ITEM: TIMESTEP
        if( !ReadLine(line,fp,iline) )
        {
            std::cerr << "EOF" << std::endl;
            return false;
        }
        
        // read timestep
        if( !ReadLine(line,fp,iline) )
        {
            throw exception("%u: missing time step",iline);
        }
        
        runtime = strconv::to<Real>(line,"Time Step");
        
        // read ITEM: NUMBER OF ATOMS
        if( !ReadLine(line,fp,iline) )
        {
            throw exception("%u: missing NUMBER OF ATOMS",iline);
        }
        
        // read #ATOMS
        if(!ReadLine(line,fp,iline))
        {
            throw exception("%u: missing #ATOMS", iline);
        }
        
        const unsigned na = strconv::to<size_t>(line,"#ATOMS");
        //std::cerr << "#ATOMS=" << na << std::endl;
        
        // skip what I don't need
        for(unsigned i=1;i<=5;++i)
        {
            if( !ReadLine(line,fp,iline) )
                throw exception("%u: Missing Header Line %u/5", iline, i);
        }
        
        // parse...
        Frame &frame = *this;
        frame.prepare(na);
        
        switch(mode)
        {
            case LoadLiquid:
                for(unsigned i=1;i<=na;++i)
                {
                    if( !ReadLine(line, fp, iline) )
                        throw exception("%u: missing liquid atom #%u", iline, i);
                    Atom atom;
                    atom.parse_full(line,iline);
                    frame.push_back(atom);
                }
                break;
                
            case LoadMembrane:
                for(unsigned i=1;i<=na;++i)
                {
                    if( !ReadLine(line, fp, iline) )
                        throw exception("%u: missing membrane atom #%u", iline, i);
                    Atom atom;
                    atom.parse_full(line,iline);
                    frame.push_back(atom);
                }
                break;

        }
        
        // a little performance counting
        const Real       t_end = chrono.query();
        const ptrdiff_t  p_end = fp.tell();
        const ptrdiff_t  nread = p_end - p_ini;
        const Real       dt    = t_end - t_ini;
        
        frame.bps = (nread/dt)*1e-6;
        
        hsort(frame, Atom::Compare );
        
        return true;
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Frame);
};


////////////////////////////////////////////////////////////////////////////////
//
// Main Program
//
////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[])
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
