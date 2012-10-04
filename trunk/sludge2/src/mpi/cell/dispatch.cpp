#include "../cell.hpp"
#include "yocto/spade/variables.hpp"


void Cell:: dispatch( const mpi &MPI )
{
    MPI.PrintfI(stderr, "layout: (%ld,%ld) -> (%ld,%ld) | y:%g -> %g | ymin=%g, ymax=%g\n", lower.x, lower.y,upper.x,upper.y,Y[lower.y],Y[upper.y],ymin,ymax);
    MPI.Printf0(stderr, "\tdispatch %u bubbles...\n", unsigned(bubbles.count()));
    bubbles.dispatch(MPI);
    
    MPI.Printf0(stderr, "\tlocate spots...\n");
    for( Bubble *b = bubbles.first(); b;b=b->next)
    {
        b->locate_spots(ymin, ymax);
    }
    
    MPI.Printf0(stderr, "\tsegmentation...\n");
    segmenter.process(bubbles);
    
    
    MPI.Printf0(stderr, "\tbuild bubble field...\n");
    segmenter.build_bubbles_field(B);
    

    
    MPI.Printf0(stderr, "\t\tsync bubble field...\n");
    save_outB( vformat("core-b%d.%d.dat",MPI.CommWorldSize,MPI.CommWorldRank));
    sync1(MPI,B);

#if 1
    segmenter.save( vformat("j%d.%d.dat",MPI.CommWorldSize,MPI.CommWorldRank));
    save_outB( vformat("sync-b%d.%d.dat",MPI.CommWorldSize,MPI.CommWorldRank));
#endif

    
    MPI.Printf0(stderr,"\tbuilding effective pressure...\n");
    segmenter.build_effective_pressure(B, P, Penter, Pleave);
    
    save_effectiveY("core");
    sync1(MPI,Penter);
    sync1(MPI,Pleave);
    save_effectiveY("sync");
    
    // no need to sync effective pressure, computed locally !!!
    
    //MPI.Printf0(stderr,"\t\tsync effective pressure...\n");
    
    
    //save_effective("eff-core.vtk");
    
    //sync1(MPI,Penter);
    //sync1(MPI,Pleave);
    
    //save_effective("eff-sync.vtk");

    
#if 0
    ios::ocstream fp("pvert-sync.dat",false);
    for( unit_t i=X.lower;i<=X.upper;++i)
    {
        fp("@i=%d\n", i);
        for( unit_t j=Y.lower;j<=Y.upper;++j)
        {
            if( Penter[j][i].y > 0 )
                fp("Penter[%d][%d].y=%g\n", j, i, Penter[j][i].y);
            if( Pleave[j][i].y > 0 )
                fp("Pleave[%d][%d].y=%g\n", j, i, Pleave[j][i].y);
        }
    }
#endif
    
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

    const Workspace &wksp = *this;
    const string     title = "effective";
    vtk.save<Layout,Real>(filename, title, wksp, pvar, outline);
    
}
