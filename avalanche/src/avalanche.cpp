#include "yocto/spade/array2d.hpp"
#include "yocto/spade/workspace.hpp"
#include "yocto/spade/rmesh.hpp"
#include "yocto/spade/vtk/writer.hpp"
#include "yocto/spade/variables.hpp"

#include "yocto/fs/vfs.hpp"
#include "yocto/fs/local-fs.hpp"

#include "yocto/exception.hpp"
#include "yocto/code/rand32.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/sort/heap.hpp"
#include "yocto/math/types.hpp"

#include <cstring>

using namespace yocto;
using namespace spade;
using namespace math;

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
        Y_SPADE_FIELD(Fields, "A", IndexArray); // 0 => water, otherwise: particle
        Y_SPADE_FIELD(Fields, "G", IndexArray); // Growing status
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};

static inline bool if_is_vtk( const vfs::entry &ep ) throw()
{
    return ep.has_extension("vtk");
}


////////////////////////////////////////////////////////////////////////////////
//
// Workspace: data + mesh
//
////////////////////////////////////////////////////////////////////////////////

#define ACTIVE_BUBBLE -1

#define AT_RIGHT  0x01
#define AT_TOP    0x02
#define AT_LEFT   0x04
#define AT_BOTTOM 0x08


class Workspace : public Parameters, public WorkspaceBase
{
public:
    Real           E[16];     //!< Energy decreasing cost
    Real           lambda;    //!< scaling factor for E0-E
    const  Real    rho;       //!< anisotropy factor
    Axis          &X;         //!< X axis
    Axis          &Y;         //!< Y axis
    IndexArray    &A;         //!< Array of state
    IndexArray    &G;         //!< Array of Growing Particles
    size_t         particles; //!< #particles
    size_t         bubbles;   //!< #bubbles
    const size_t   Nx;
    const size_t   Ny;
    const size_t   Nc;        //!< (Nx-2) * (Ny-2)
    const Real     frequency; //!< one spontaneous bubble every 1/frequency step per slot
    const Real     E0;
    uniform_generator<double>  ran;
    vector<size_t> Size;  //!< Size of each particle
    vtk_writer     vtk;
    variables      var;
    string         outdir;
    vfs           &fs;
    
    explicit Workspace(const Layout &l,
                       const double Rho = 0.5,
                       const double Gain = 1
                       ) :
    WorkspaceBase(l,Fields,Ghosts),
    E(),
    lambda(2*Log(max_of<Real>(1,Gain))),
    rho(clamp<Real>(0,Rho,1)),
    X( mesh.X() ),
    Y( mesh.Y() ),
    A( (*this)["A"].as<IndexArray>() ),
    G( (*this)["G"].as<IndexArray>() ),
    particles(0),
    bubbles(0),
    Nx( l.width.x ),
    Ny( l.width.y ),
    Nc( (Nx-2) * (Ny-2) ),
    frequency(1e-4),
    E0( -Log(frequency) ),
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
        
        fs.as_directory(outdir);
        fs.create_dir(outdir,true);
        fs.remove_files(outdir, if_is_vtk);
        
        memset(E,0,sizeof(E));
        
        const Real scale2 = 1.0/2;
        const Real scale3 = 1.0/3;
        const Real scale4 = 1.0/4;
        
        //----------------------------------------------------------------------
        // Zero Neighbors: 1 case
        //----------------------------------------------------------------------
        E[0] = 0; //!< spontaneous
        
        //----------------------------------------------------------------------
        // One Neighbor: 4 cases
        //----------------------------------------------------------------------
        const Real h_rho = rho;
        const Real v_rho = 1-h_rho;
        E[AT_LEFT]   = h_rho;
        E[AT_RIGHT]  = h_rho;
        E[AT_TOP]    = v_rho;
        E[AT_BOTTOM] = v_rho;
        
        //----------------------------------------------------------------------
        // Two Neighbors: 6 cases
        //----------------------------------------------------------------------
        E[AT_LEFT|AT_RIGHT]   = scale2 * ( E[AT_LEFT]  + E[AT_RIGHT]  );
        E[AT_LEFT|AT_TOP]     = scale2 * ( E[AT_LEFT]  + E[AT_TOP]    );
        E[AT_LEFT|AT_BOTTOM]  = scale2 * ( E[AT_LEFT]  + E[AT_BOTTOM] );
        E[AT_RIGHT|AT_TOP]    = scale2 * ( E[AT_RIGHT] + E[AT_TOP]    );
        E[AT_RIGHT|AT_BOTTOM] = scale2 * ( E[AT_RIGHT] + E[AT_BOTTOM] );
        E[AT_TOP|AT_BOTTOM]   = scale2 * ( E[AT_TOP]   + E[AT_BOTTOM] );
        
        //----------------------------------------------------------------------
        // Three Neighbors: 4 cases case
        //----------------------------------------------------------------------
        E[AT_LEFT|AT_RIGHT|AT_TOP]    = scale3 * ( E[AT_LEFT] + E[AT_RIGHT]  + E[AT_TOP]    );
        E[AT_LEFT|AT_RIGHT|AT_BOTTOM] = scale3 * ( E[AT_LEFT] + E[AT_RIGHT]  + E[AT_BOTTOM] );
        E[AT_TOP|AT_BOTTOM|AT_LEFT]   = scale3 * ( E[AT_TOP]  + E[AT_BOTTOM] + E[AT_LEFT]   );
        E[AT_TOP|AT_BOTTOM|AT_RIGHT]  = scale3 * ( E[AT_TOP]  + E[AT_BOTTOM] + E[AT_RIGHT]  );
        
        
        //----------------------------------------------------------------------
        // Four Neighbors: one case
        //----------------------------------------------------------------------
        E[AT_LEFT|AT_RIGHT|AT_BOTTOM|AT_TOP] = scale4*(E[AT_TOP]  + E[AT_BOTTOM] + E[AT_RIGHT] + E[AT_LEFT] );
        
        
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
        A.ldz();
        G.ldz();
        Size.free();
        
    CREATE_BUBBLES:
        //-- let us make initial_bubbles
        for(size_t j=2;j<Ny;++j)
        {
            for(size_t i=2;i<Nx;++i)
            {
                if( A[j][i]<=0 && ran() <= level )
                {
                    ++particles;
                    ++bubbles;                 // initially: one particle = one bubble
                    A[j][i] = particles;       // a new particle
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
    
    
    
    void check_owner( size_t i, size_t j, size_t *owner, size_t &owners, size_t &weight) const throw()
    {
        assert(owner);
        assert(owners<4);
        
        const size_t particle_index = A[j][i];
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
                if(A[j][i]>0)
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
                assert(flag<16);
                
                //--------------------------------------------------------------
                // initialize cost
                //--------------------------------------------------------------
                const Real alpha = ran();
                const Real dE    = E0 - lambda * E[flag];
                
                
                if(alpha < Exp(-dE) )
                {
                    if(flag>0)
                    {
                        //------------------------------------------------------
                        // find the particle
                        //------------------------------------------------------
                        assert(owners>0);
                        const size_t p = owner[ ran.lt(owners) ];
                        assert(p>=1);
                        assert(p<=particles);
                        assert(Size.size()==particles);
                        A[j][i] = p;
                        G[j][i] = p;
                        ++Size[p];
                    }
                    else
                    {
                        //------------------------------------------------------
                        // create a particle
                        //------------------------------------------------------
                        ++particles;
                        ++bubbles;
                        Size.push_back(1);
                        assert(Size.size()==particles);
                        G[j][i] = particles;
                        A[j][i] = particles;
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
        size_t        Nx = 200;
        size_t        Ny = 200;
        const  Layout LL( Coord(1,1), Coord(Nx,Ny) );
        Workspace     W(LL,0.45,3500);
        const size_t  ini      = 0;
        const size_t  iter_max = 400;
        
        std::cerr << "W.frequency = " << W.frequency << std::endl;
        std::cerr << "W.lambda    = " << W.lambda    << std::endl;
        std::cerr << "W.E0        = " << W.E0        << std::endl;
        std::cerr << "W.Ncore     = " << W.Nc        << std::endl;
        
        W.reset(ini);
        W.save(0);
        for(size_t iter=1;iter<=iter_max;++iter)
        {
            W.step();
            W.save(iter);
        }
        
        std::cerr << "#particles=" << W.particles << std::endl;

        hsort(W.Size);
        W.Size.reverse();
        {
            ios::ocstream fp("dist.dat",false);
            for(size_t i=1; i<= W.Size.size(); ++i )
            {
                fp("%g %g\n", double(i), double(W.Size[i]));
            }
        }
        return 0;
    }
    catch(...)
    {
        std::cerr << "unhandled exception in " << prog << std::endl;
    }
    return 1;
}
