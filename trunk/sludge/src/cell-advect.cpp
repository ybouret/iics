#include "cell.hpp"


void Cell:: advect(Real dt)
{
    for( Bubble *bubble = bubbles.first(); bubble;bubble=bubble->next)
    {
        for( Spot *spot = bubble->spots.head;spot;spot=spot->next)
        {
            Tracer *p  = spot->handle;
            p->vertex += dt * spot->U;
        }
    }
}