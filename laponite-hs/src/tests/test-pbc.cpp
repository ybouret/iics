#include "../types.hpp"

static
void show( Real y, const PBC &pbc )
{
    
    std::cerr << y << " [" << pbc.L << "]=" << pbc(y) << std::endl;
}

int main( int argc, char *argv[] )
{
    
    const Real L=3.3;
    const PBC  pbc(L);
    
    show(L, pbc);
    show(-L,pbc);
    show( 3.5, pbc );
    show(-3.5,pbc);
    
    show(L/2,pbc);
    show(-L/2,pbc);
    
    return 0;
}