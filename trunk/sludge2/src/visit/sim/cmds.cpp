#include "../simulation.hpp"


void Simulation:: perform(const string &cmd)
{
    if( cmd == "raz" )
    {
        initialize();
        return;
    }
}

