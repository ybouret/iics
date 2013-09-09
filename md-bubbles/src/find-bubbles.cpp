#include "yocto/fs/local-fs.hpp"
#include "yocto/math/v3d.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/string/tokenizer.hpp"
#include "yocto/exception.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/math/alg/delaunay.hpp"
#include "yocto/auto-ptr.hpp"
#include "yocto/sort/heap.hpp"
#include "yocto/sys/wtime.hpp"
#include <iostream>

using namespace yocto;
using namespace math;

typedef double    Real;
typedef v3d<Real> V3D;

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
    
    inline void parse( const string &line, unsigned iline)
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

typedef vector<Atom>      Atoms;
typedef vector<V3D>       Vertices;
typedef vector<iTriangle> Triangles;

class Frame : public Atoms
{
public:
    Frame *next;
    Frame *prev;
    Real   v_ave;
    Real   v_sig;
    Real   bps;
    Real   fps;
    wtime  chrono;
    
    explicit Frame() throw() :
    Atoms(),
    next(0),
    prev(0),
    v_ave(0),
    v_sig(0),
    bps(0),
    fps(0)
    {
        chrono.start();
    }
    
    void prepare(size_t na)
    {
        free();
        ensure(na);
        v_ave = 0;
        v_sig = 0;
        bps   = 0;
    }
    
    void find_gas(Vertices &gas, Real vmin ) const
    {
        gas.free();
        for(size_t j=size();j>0;--j)
        {
            const Atom &atom = (*this)[j];
            if(atom.voronoi>=vmin)
            {
                gas.push_back( atom.r );
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
    
    
    bool load_next(ios::icstream &fp, unsigned &iline)
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
        
        Real vor_max = 0;
        for(unsigned i=1;i<=na;++i)
        {
            if( !ReadLine(line, fp, iline) )
                throw exception("%u: missing atom #%u", iline, i);
            Atom atom;
            atom.parse(line,iline);
            frame.push_back(atom);
            const Real vor = atom.voronoi;
            frame.v_ave += vor;
            if(vor>vor_max)
                vor_max = vor;
        }
        const Real       t_end = chrono.query();
        const ptrdiff_t  p_end = fp.tell();
        const ptrdiff_t  nread = p_end - p_ini;
        const Real dt    = t_end - t_ini;
        frame.fps    = 1.0 / dt;
        frame.bps    = (nread/dt)*1e-6;
        frame.v_ave /= na;
        //std::cerr << "vor_max=" << vor_max << std::endl;
        
        //hsort(*frame, Atom::Compare);
        
        return true;
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Frame);
};


int main( int argc, char *argv[] )
{
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        //======================================================================
        //
        // Open File
        //
        //======================================================================
        if(argc<=2)
            throw exception("usage: %s input_file vmin [#frame_max]", prog);
        
        const string  input_name = argv[1];
        const Real    vmin       = strconv::to<Real>( argv[2], "vmin");
        size_t        frame_max  = 0;
        if(argc>3)
            frame_max = strconv::to<size_t>(argv[3],"#frame_max");
        
        ios::icstream fp( input_name  );
        
        unsigned        iline = 1;
        size_t          num_frame = 0;
        
        if(frame_max>0) std::cerr << "Reading at most " << frame_max << " frames" << std::endl;
        std::cerr << "In gas <=> voronoi >= " << vmin << std::endl;
        Frame    frame;
        Vertices gas;
        while( frame.load_next(fp,iline) )
        {
            ++num_frame;
            
            //------------------------------------------------------------------
            // process
            //------------------------------------------------------------------
            frame.find_gas(gas, vmin);
            std::cerr << "#" << num_frame
            << " <voronoi>="
            << frame.v_ave
            << " : #gas=" << gas.size()
            << " fps=" << frame.fps
            << " @" << frame.bps << " MB/s"
            << std::endl;
            
            if(frame_max>0 && num_frame>=frame_max)
                break;
        }
        std::cerr << std::endl;
        
        return 0;
    }
    catch( const exception &e )
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "Unhandled Exception in " << prog << std::endl;
    }
    return 1;
}
