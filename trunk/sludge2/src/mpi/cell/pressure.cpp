#include "../cell.hpp"

void Cell:: init_pressure(const mpi &MPI )
{
    //! bubble
    for( const Marker *m = segmenter.markers.head;m;m=m->next)
    {
        P[m->inside.y][m->inside.x] = m->bubble->pressure;
    }
    
    //! boundary
    
    //! sync pressure
    sync1(MPI,P);
    
}