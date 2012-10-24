#include "../cell.hpp"
#include "yocto/spade/variables.hpp"


void Cell:: dispatch( const mpi &MPI )
{
    //MPI.PrintfI(stderr, "layout: (%ld,%ld) -> (%ld,%ld) | y:%g -> %g | ymin=%g, ymax=%g\n", lower.x, lower.y,upper.x,upper.y,Y[lower.y],Y[upper.y],ymin,ymax);
    
    //--------------------------------------------------------------------------
    // dispatch bubbles
    //--------------------------------------------------------------------------
    MPI.Printf0(stderr, "\tdispatch %u bubbles...\n", unsigned(bubbles.count()));
    bubbles.dispatch(MPI);
    
    //--------------------------------------------------------------------------
    // locate who is in our domain
    //--------------------------------------------------------------------------
    MPI.Printf0(stderr, "\tlocate spots...\n");
    for( Bubble *b = bubbles.first(); b;b=b->next)
    {
        b->locate_spots(ymin, ymax);
    }
    
    
    //--------------------------------------------------------------------------
    // find the interior of bubbles
    //--------------------------------------------------------------------------
    MPI.Printf0(stderr, "\tsegmentation...\n");
    segmenter.process(bubbles);
    //segmenter.save( vformat("core-j%d.%d.dat",MPI.CommWorldSize,MPI.CommWorldRank));
    
    //--------------------------------------------------------------------------
    // build the local bubble field
    //--------------------------------------------------------------------------
    MPI.Printf0(stderr, "\tbuild bubble field...\n");
    segmenter.build_bubbles_field(B);
    
    //--------------------------------------------------------------------------
    // synchronize the bubble field
    //--------------------------------------------------------------------------
    MPI.Printf0(stderr, "\tsync bubble field...\n");
    sync1(MPI,B);
    
    save_outB("b.dat");
    
    MPI.Printf0(stderr, "\tbuild bulk field...\n");
    build_bulk();
    
    //--------------------------------------------------------------------------
    // Junctions PBC
    //--------------------------------------------------------------------------
    MPI.Printf0(stderr, "\tdispatch vertical junctions...\n");
    segmenter.dispatch_vertical_junctions(MPI, *this);
    
    //--------------------------------------------------------------------------
    // Effective Pressure
    //--------------------------------------------------------------------------
    MPI.Printf0(stderr,"\tbuilding effective pressure...\n");
    segmenter.build_effective_pressure(B, P, Penter, Pleave, bubbles.gamma);
    
    
}


static inline void save_field_Y( const string &filename,
                                const VertexArray &A,
                                const Array1D     &X,
                                const Array1D     &Y)
{
    ios::ocstream fp(filename,false);
    for( unit_t j=A.lower.y;j<=A.upper.y;++j)
    {
        for( unit_t i=A.lower.x;i<=A.upper.x;++i)
        {
            const Real a = A[j][i].y;
            if(a>0)
            {
                fp("%g %g\n", X[i],Y[j]);
            }
        }
    }
}

void Cell:: save_effectiveY(const string &prefix) const
{
    {
        const string filename = prefix + "-enter.dat";
        save_field_Y(filename,Penter,X,Y);
    }
    
    {
        const string filename = prefix + "-leave.dat";
        save_field_Y(filename,Pleave,X,Y);
    }
}

void Cell:: save_effective( const string &filename) const
{
    
    variables pvar;
    pvar.append("Penter");
    pvar.append("Pleave");
    
    const Workspace &wksp  = *this;
    const string     title = "effective";
    vtk.save<Layout,Real>(filename, title, wksp, pvar, outline);
    
}
