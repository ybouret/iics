#include "bubbles.hpp"




void Bubbles:: check_geometries_within( Real y_lo, Real y_hi )
{
    for( Bubble *bubble = first(); bubble; bubble=bubble->next )
    {
        bubble->collect_spots_within(y_lo, y_hi);
        if( bubble->active )
        {
            bubble->compute_geometry();
        }
    }
}
