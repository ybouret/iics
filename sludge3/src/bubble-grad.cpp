#include "bubble.hpp"

void Bubble:: bracket( const Marker *m, const Junction **Jprev, const Junction **Jnext) const
{
    assert(m);
    assert(markers.owns(m));
    assert(Jprev);
    assert(Jnext);
    
}

void Bubble:: normal_grad( const Marker *m, const Array &B, const Array &P)
{
    const Junction *Jprev=0, *Jnext=0;
    bracket(m, &Jprev, &Jnext);

}
