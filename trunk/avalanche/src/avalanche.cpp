#include "yocto/spade/array2d.hpp"
#include "yocto/spade/workspace.hpp"
#include "yocto/spade/rmesh.hpp"
#include "yocto/spade/vtk/writer.hpp"
#include "yocto/spade/variables.hpp"

#include "yocto/fs/vfs.hpp"
#include "yocto/fs/local-fs.hpp"

#include "yocto/exception.hpp"
#include "yocto/code/rand32.hpp"

using namespace yocto;
using namespace spade;

////////////////////////////////////////////////////////////////////////////////
//
// Defining types
//
////////////////////////////////////////////////////////////////////////////////

typedef coord2D                        Coord;
typedef double                         Real;
typedef layout2D                       Layout;
typedef array1D<Real>                  Axis;
typedef array2D<Real>                  Array;
typedef array2D<unit_t>                IndexArray;
typedef fields_setup<Layout>           FieldsSetup;
typedef ghosts_setup                   GhostsSetup;
typedef workspace<layout2D,rmesh,Real> WorkspaceBase;


class URand : public rand32_kiss
{
public:
    explicit URand() throw()
    {
        wseed();
    }
    
    virtual ~URand() throw()
    {
    }
    
    inline Real operator()(void) throw()
    {
        return get<Real>();
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(URand);
};


////////////////////////////////////////////////////////////////////////////////
//
// Parameters for the Workspace
//
////////////////////////////////////////////////////////////////////////////////

class Parameters
{
public:
    FieldsSetup Fields;
    GhostsSetup Ghosts;
    explicit Parameters() :
    Fields(0),
    Ghosts()
    {
        Y_SPADE_FIELD(Fields, "A", Array);      // Array
        Y_SPADE_FIELD(Fields, "B", IndexArray); // Belonging to
        Y_SPADE_FIELD(Fields, "G", IndexArray); // Growing status
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};

////////////////////////////////////////////////////////////////////////////////
//
// Workspace: data + mesh
//
////////////////////////////////////////////////////////////////////////////////

#define ACTIVE_BUBBLE -1

static bool if_is_vtk( const vfs::entry &ep ) throw()
{
    const char *ext = ep.extension;
    return ext && strcmp(ext,"vtk") == 0;
}

class Workspace : public Parameters, public WorkspaceBase
{
public:
    Axis          &X;         //!< X axis
    Axis          &Y;         //!< Y axis
    Array         &A;         //!< Array of state
    IndexArray    &B;         //!< Array of bubbles segementation
    IndexArray    &G;         //!< Array of Growing Particles
    size_t         particles; //!< #particles
    size_t         bubbles;   //!< #bubbles
    const size_t   Nx;
    const size_t   Ny;
    const size_t   Nc;    //!< (Nx-2) * (Ny-2)
    const size_t   M;
    const Real     spontaneous_level;
    URand          ran;
    vector<size_t> Size;  //!< Size of each particle
    vtk_writer     vtk;
    variables      var;
    string         outdir;
    vfs           &fs;
    
    explicit Workspace( const Layout &l ) :
    WorkspaceBase(l,Fields,Ghosts),
    X( mesh.X() ),
    Y( mesh.Y() ),
    A( (*this)["A"].as<Array>() ),
    B( (*this)["B"].as<IndexArray>() ),
    G( (*this)["G"].as<IndexArray>() ),
    particles(0),
    bubbles(0),
    Nx( l.width.x ),
    Ny( l.width.y ),
    Nc( (Nx-2) * (Ny-2) ),
    M(2),
    spontaneous_level(1.0/(M*Nc)),
    ran(),
    Size(Nc,as_capacity),
    vtk(),
    var(),
    outdir("out/"),
    fs( local_fs::instance() )
    {
        for(size_t i=1; i <= Nx; ++i ) X[i] = i;
        for(size_t j=1; j <= Ny; ++j ) Y[j] = j;
        
        var.append("A");
        var.append("B");
        
        fs.as_directory(outdir);
        fs.create_dir(outdir,true);
        fs.remove_files(outdir, if_is_vtk);
    }
    
    void save( size_t idx ) const
    {
        const string filename = outdir + vformat("A%u.vtk", unsigned(idx));
        vtk.save(filename, "fields", *this, var, as_layout());
    }
    
    
    void reset(size_t initial_bubbles) throw()
    {
        if(initial_bubbles>Nc) initial_bubbles = Nc;
        const Real level = initial_bubbles / double(Nc);
        
        particles = 0;
        bubbles   = 0;
        //-- fill all with water
        A.ld(1);
        B.ldz();
        G.ldz();
        Size.free();
        
    CREATE_BUBBLES:
        //-- let us make initial_bubbles
        for(size_t j=2;j<Ny;++j)
        {
            for(size_t i=2;i<Nx;++i)
            {
                if( B[j][i]<=0 && ran() <= level )
                {
                    A[j][i] = 0;
                    ++particles;
                    ++bubbles;                 // initially: one particle = one bubble
                    B[j][i] = particles;       // a new particle
                    G[j][i] = ACTIVE_BUBBLE;   // is growing
                    Size.push_back(1);         // size of particle #bubble
                    if(particles>=initial_bubbles)
                        goto DONE_WITH_BUBBLES;
                }
            }
        }
        if(particles<initial_bubbles)
            goto CREATE_BUBBLES;
    DONE_WITH_BUBBLES:
        ;
    }
    
#define AT_LEFT   0x01
#define AT_RIGHT  0x02
#define AT_TOP    0x04
#define AT_BOTTOM 0x08
    
    void check_owner( size_t i, size_t j, size_t *owner, size_t &owners, size_t &weight) const throw()
    {
        assert(owner);
        assert(owners<4);
        
        const size_t particle_index = B[j][i];
        assert(particle_index>0);
        assert(particle_index<=particles);
        assert(Size.size() == particles);
        const size_t w  = Size[particle_index];
        if(w<weight)
            return; // don't change
        if(w>weight)
        {
            // take precendence
            owners = 1;
            owner[0] = particle_index;
            weight = w;
            return;
        }
        assert(w==weight);
        assert(owners>0);
        // multiple choice
        owner[owners++] = particle_index;
    }
    
    void step()
    {
        //-- first pass: find the active bubbles
        for(size_t j=2;j<Ny;++j)
        {
            const size_t jm = j-1;
            const size_t jp = j+1;
            for(size_t i=2;i<Nx;++i)
            {
                //-- a bubble remains a bubble
                if(A[j][i]<=0)
                    continue;
                
                //-- we are in the water
                const size_t im=i-1;
                const size_t ip=i+1;
                unsigned flag     = 0;
                size_t   owner[4] = { 0 };
                size_t   owners   = 0;
                size_t   weight   = 0;
                
                //--------------------------------------------------------------
                // detect configuration
                //--------------------------------------------------------------
                if( G[j][i-1] == ACTIVE_BUBBLE )
                {
                    flag |= AT_LEFT;
                    check_owner(im, j, owner, owners, weight);
                }
                
                if( G[j][ip] == ACTIVE_BUBBLE )
                {
                    flag |= AT_RIGHT;
                    check_owner(ip, j, owner, owners, weight);
                }
                
                if( G[jm][i] == ACTIVE_BUBBLE )
                {
                    flag |= AT_BOTTOM;
                    check_owner(i, jm, owner, owners, weight);
                }
                
                if( G[jp][i] == ACTIVE_BUBBLE )
                {
                    flag |= AT_TOP;
                    check_owner(i, jp, owner, owners, weight);
                }
                
                //--------------------------------------------------------------
                // initialize cost
                //--------------------------------------------------------------
                const Real r = ran();
                if(flag>0)
                {
                    
                }
                else
                {
                    if(r<spontaneous_level)
                    {
                        std::cerr << "New Particle" << std::endl;
                        ++particles;
                        ++bubbles;
                        Size.push_back(1);
                        B[j][i] = particles;
                        G[j][i] = particles;
                        A[j][i] = 0;
                    }
                }
                
                
            }
        }
        
        //-- second pass: regularize status
        for(size_t j=2;j<Ny;++j)
        {
            for(size_t i=2;i<Nx;++i)
            {
                // turn off active bubble
                if( G[j][i] == ACTIVE_BUBBLE)
                {
                    G[j][i] = 0;
                    continue;
                }
                
                // turn on newly created bubbles
                if(G[j][i]>0)
                {
                    G[j][i] = ACTIVE_BUBBLE;
                }
            }
            
        }
    }
    
    virtual ~Workspace() throw()
    {
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Workspace);
};

////////////////////////////////////////////////////////////////////////////////
//
// main
//
////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        size_t        Nx = 100;
        size_t        Ny = 80;
        const  Layout LL( Coord(1,1), Coord(Nx,Ny) );
        Workspace     W(LL);
        const size_t  ini = 10;
        const size_t  iter_max = 100;
        
        W.reset(ini);
        
        W.save(0);
        
        for(size_t iter=1;iter<=iter_max;++iter)
        {
            W.step();
            W.save(iter);
        }
        
        std::cerr << "#particles=" << W.particles << "/" << ini + double(iter_max)/W.M << std::endl;
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "unhandled exception in " << prog << std::endl;
    }
    return 1;
}
