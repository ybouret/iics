#include "cell.hpp"


void Cell:: assemble_all()
{
    bubbles.assemble_all(MPI);
}
