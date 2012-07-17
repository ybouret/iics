#include "../cell.hpp"

#include "yocto/ios/ocstream.hpp"


static inline void save_bubble( const Bubble &b )
{
    static int id = 0;
    ios::ocstream fp( vformat("bubble%d.curve",id++), false );
    const Point *p = b.root;
    fp("#curve\n");
    for(size_t i=b.size;i>0;--i,p=p->next )
    {
        fp("%g %g\n", p->vertex.x, p->vertex.y);
    }
    fp("%g %g\n", p->vertex.x, p->vertex.y);
    fp("\n");
}

int main( int argc, char *argv[] )
{
    try 
    {
        AleaInit();
        mpi &MPI = mpi::init( &argc, &argv);
        
        Cell cell(10,20,5.0,6.0,MPI);
        MPI.Printf(stderr, "rank %d> SubRegion y: %g -> %g\n", MPI.CommWorldRank, cell.SubRegion.vmin.y, cell.SubRegion.vmax.y );
        
        
        if( MPI.IsMaster )
        {
            Bubble *b = cell.bubbles.create();
            b->map_circle(V2D(0,0), 1.2);
        }
        
        
        
        if( MPI.IsMaster )
        {
            cell.master_update_topologies();
            save_bubble( *cell.bubbles.first() );
        }
        
        
        for( size_t iter=0; iter<100;++iter )
        {
            
            
            //! rank 0 -> everybody
            cell.dispatch_bubbles(MPI);
            
            //! move bubble
            for( Bubble *b = cell.bubbles.first(); b; b=b->next )
            {
                for( Spot *spot = b->spots.head; spot != NULL; spot = spot->next )
                {
                    spot->point->vertex.y += 0.01;
                    spot->point->vertex.x *= 1.05;
                }
            }
            
            //! and everybody -> rank 0
            cell.assemble_bubbles(MPI);
            MPI.Printf0( stderr, "." );
            if( MPI.IsMaster )
            {
                cell.master_update_topologies();
                save_bubble( *cell.bubbles.first() );
            }
        }
        MPI.Printf0(stderr, "\n");
        
        
        
        return 0;
    }
    catch(...)
    {
        std::cerr << "Unhandled exception" << std::endl;
    }
    return 1;
}
