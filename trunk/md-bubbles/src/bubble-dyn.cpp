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
#include "yocto/memory/buffers.hpp"

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
    
    // read membrane atom
    inline void parse_fast( const string &line, unsigned iline)
    {
        tokenizer tkn(line);
        
        //======================================================================
        // read id
        //======================================================================
        id = Next<size_t>(tkn,iline,"id");
        (void)Next<size_t>(tkn,iline,"type");
        
        //======================================================================
        // read position
        //======================================================================
        r.x = Next<Real>(tkn, iline, "x");
        r.y = Next<Real>(tkn, iline, "y");
        r.z = Next<Real>(tkn, iline, "z");
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
    Vertices       gas;    //(1024,as_capacity);
    Sorted         gid;    //(1024,as_capacity);
    Triangles      trlist; //(1024,as_capacity);
    vector<size_t> hull;   //(1024,as_capacity); // for hull
    
    explicit Frame() throw() :
    Atoms(),
    next(0),
    prev(0),
    runtime(0),
    bps(0),
    chrono(),
    gas(),
    gid(),
    trlist(),
    hull()
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
    
    inline void find_gas(Real vmin)
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
                    atom.parse_fast(line,iline);
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
    
    
    void process_gas(const double Vmax,
                     const string &output_name,
                     const string &xyz_name)
    {
        //-- find the gas particles
        find_gas(Vmax);
        delaunay<Real>::build(trlist, gas);
        delaunay_hull(hull, trlist);
        
        //-- compute and save the area
        const Real area = delaunay<Real>::area(hull,gas);
        
        {
            ios::ocstream output(output_name,true);
            output("%g %g %g\n",  runtime, area, double(gas.size()));
        }
        
        //-- is there is some gas=>xyz
        const size_t ng = gas.size();
        if(ng>0)
        {
            ios::ocstream xyz(xyz_name,true);
            xyz("%u\n", unsigned( gas.size() ) );
            xyz("t=%g\n",runtime);
            for(size_t i=1;i<=gas.size();++i)
            {
                const char *name = "H";
                //if( gid.search(i) ) name = "He";
                //const Atom &atom = frame[i];
                const Atom &atom = (*this)[ gid[i] ];
                xyz("%s %.4e %.4e %.4e\n", name, atom.r.x, atom.r.y, atom.r.z);
            }

        }
        
    }
    
    
    void process_membrane(const string &output_name,
                          const string &xyz_name)
    {
        //-- save as XYZ
        {
            ios::ocstream xyz(xyz_name,true);
            xyz("%u\n", unsigned( size() ) );
            xyz("t=%g\n",runtime);
            for(size_t i=1;i<=size();++i)
            {
                const char *name = "He";
                const Atom &atom = (*this)[ i ];
                xyz("%s %.4e %.4e %.4e\n", name, atom.r.x, atom.r.y, atom.r.z);
            }
        }
        
        //-- compute barycenter and save it...
        {
            ios::ocstream out(output_name,true);
            V3D G;
            for(size_t i=size();i>0;--i)
            {
                G += (*this)[i].r;
            }
            G *= 1/Real(size());
            out("%g %.4e %.4e %.4e\n", runtime, G.x, G.y, G.z);
        }
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
        vfs &fs = local_fs::instance();
        Real Vmax = 0;
        
        //======================================================================
        //
        // Parsing arguments
        //
        //======================================================================
        const string  input_name   = argv[1];
        LoadMode      load_mode    = LoadLiquid;
        string        output_name  = input_name;
        string        xyz_name     = input_name;
        vfs::change_extension(xyz_name,"xyz");
        switch(argc)
        {
            case 2:
                load_mode = LoadMembrane;
                std::cerr << "-- <Loading Membrane from " << vfs::get_base_name(input_name) << ">" << std::endl;
                fs.change_extension(output_name, "memb.dat");
                break;
                
            case 3:
                load_mode = LoadLiquid;
                std::cerr << "-- <Loading Liquid from .../" << vfs::get_base_name(input_name) << ">" << std::endl;
                Vmax = strconv::to<Real>(argv[2],"VoronoiCutOff");
                fs.change_extension(output_name, "area.dat");
                break;
                
            default:
                throw exception("usage: %s data [ VoronoiCutOff ]",prog);
        }
        std::cerr << "-- <Saving data into '" << vfs::get_base_name(output_name) << "'>" << std::endl;
        std::cerr << "-- <Saving XYZ  into '" << vfs::get_base_name(xyz_name) << "'>" << std::endl;
        
        ios::ocstream::overwrite(output_name);
        ios::ocstream::overwrite(xyz_name);
        
        //======================================================================
        //
        // Preparing I/O
        //
        //======================================================================
        ios::icstream                             fp( input_name  );
        memory::buffer_of<uint8_t,memory::global> iobuf(32*1024*1024);
        fp.bufferize(iobuf);
        
        
        //======================================================================
        //
        // Processing
        //
        //======================================================================
        Frame        frame;
        unsigned     iline       = 1;
        unsigned     num_frame   = 0;
        double       average_bps = 0;
        const size_t every       = 10;
        
        std::cerr.flush();
        while( frame.load_next(fp,iline,load_mode) )
        {
            //------------------------------------------------------------------
            // got a frame
            //------------------------------------------------------------------
            ++num_frame;
            average_bps += frame.bps;
            if(0==(num_frame%every))
            {
                average_bps /= every;
                fprintf(stderr,"-- Read @ %12.6f MBytes/s\n", average_bps);
                average_bps = 0;
            }
            
            //------------------------------------------------------------------
            // process according to file type
            //------------------------------------------------------------------
            switch (load_mode)
            {
                case LoadLiquid:
                    frame.process_gas(Vmax,output_name,xyz_name);
                    break;
                    
                case LoadMembrane:
                    frame.process_membrane(output_name, xyz_name);
                    break;
            }
        }
        std::cerr.flush();
        std::cerr << std::endl;
        
        return 0;
    }
    catch(const exception &e)
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
