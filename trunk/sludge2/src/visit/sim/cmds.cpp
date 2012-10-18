#include "../simulation.hpp"


void Simulation:: perform(const string &cmd, const array<string> &args)
{
    if( cmd == "raz" )
    {
        initialize();
        return;
    }
}

