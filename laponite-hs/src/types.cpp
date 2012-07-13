#include "types.hpp"


#include "yocto/code/rand.hpp"

void AleaInit() throw()
{
    _rand.wseed();
}

Real Alea() throw()
{
    return alea<Real>();
}