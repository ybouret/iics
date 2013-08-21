#include "yocto/spade/array2d.hpp"
#include "yocto/spade/workspace.hpp"
#include "yocto/spade/rmesh.hpp"
#include "yocto/spade/vtk/writer.hpp"
#include "yocto/spade/variables.hpp"

#include "yocto/fs/vfs.hpp"
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
        Y_SPADE_FIELD(Fields, "A", Array); // Array
        Y_SPADE_FIELD(Fields, "B", IndexArray); // Belonging to
        Y_SPADE_FIELD(Fields, "G", IndexArray); // Growing
    }
    
private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Parameters);
};

////////////////////////////////////////////////////////////////////////////////
//
// Workspace: data + mesh
//
////////////////////////////////////////////////////////////////////////////////

#define PASSIVE_BUBBLE 0
#define ACTIVE_BUBBLE -1

class Workspace : public Parameters, public WorkspaceBase
{
public:
    Axis        &X;       //!< X axis
    Axis        &Y;       //!< Y axis
    Array       &A;       //!< Array of state
    IndexArray  &B;       //!< Array of bubbles segementation
    IndexArray  &G;       //!< Array of Growing Particles
    size_t       bubbles; //!< #bubbles
    Real         h_energy;
    Real         v_energy;
    const size_t Nx;
    const size_t Ny;
    const size_t Nc;   //!< (Nx-2) * (Ny-2)
    URand        ran;  
    vtk_writer   vtk;
    
    explicit Workspace( const Layout &l ) :
    WorkspaceBase(l,Fields,Ghosts),
    X( mesh.X() ),
    Y( mesh.Y() ),
    A( (*this)["A"].as<Array>() ),
    B( (*this)["B"].as<IndexArray>() ),
    G( (*this)["G"].as<IndexArray>() ),
    bubbles(0),
    h_energy(1),
    v_energy(1),
    Nx( l.width.x ),
    Ny( l.width.y ),
    Nc( (Nx-2) * (Ny-2) )
    {
        for(size_t i=1; i <= Nx; ++i ) X[i] = i;
        for(size_t j=1; j <= Ny; ++j ) Y[j] = j;
    }
    
    
    
    void reset(size_t initial_bubbles) throw()
    {
        if(initial_bubbles>Nc) initial_bubbles = Nc;
        const Real level = initial_bubbles / double(Nc);

        bubbles = 0;
        //-- fill all with water
        A.ld(1);
        B.ldz();
        G.ldz();
        
    CREATE_BUBBLES:        
        //-- let us make initial_bubbles
        for(size_t j=2;j<Ny;++j)
        {
            for(size_t i=2;i<Nx;++i)
            {
                if( B[j][i]<=0 && ran() <= level )
                {
                    A[j][i] = 0;
                    ++bubbles;
                    B[j][i] = bubbles;       // a new particle
                    G[j][i] = ACTIVE_BUBBLE; // is growing
                    if(bubbles>=initial_bubbles)
                        goto DONE_WITH_BUBBLES;
                }
            }
        }
        if(bubbles<initial_bubbles)
            goto CREATE_BUBBLES;
    DONE_WITH_BUBBLES:
        ;
    }
    

    void step()
    {
        //-- first pass: find the active bubbles
        for(size_t j=2;j<Ny;++j)
        {
            for(size_t i=2;i<Nx;++i)
            {
                Real E = 0;
                //==============================================================
                // write a decent activation model !
                //==============================================================
                if(G[j][i-1] == ACTIVE_BUBBLE || G[j][i+1] == ACTIVE_BUBBLE)
                    E += h_energy;
               
                if(G[j-1][i] == ACTIVE_BUBBLE || G[j+1][i] == ACTIVE_BUBBLE)
                    E += v_energy;
                
                //==============================================================
                // change
                //==============================================================
                
                
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
        
        variables var;
        var.append("A");
        var.append("B");
        
        W.reset(10);
        
        W.vtk.save("ini.vtk", "init", W, var, W.as_layout());
        
        W.step();
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "unhandled exception in " << prog << std::endl;
    }
    return 1;
}
