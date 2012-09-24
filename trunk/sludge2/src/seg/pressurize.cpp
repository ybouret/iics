#include "../segmenter.hpp"

void Segmenter:: pressurize(Array &P) const
{

    for( const Marker *m = markers.head;m;m=m->next)
    {
        const Bubble *bubble = m->bubble; assert(bubble);
        P[m->inside.y][m->inside.x] = bubble->pressure;
    }

}