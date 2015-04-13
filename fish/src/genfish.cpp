#include "fish.hpp"
#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"

YOCTO_PROGRAM_START()
{
    Lua::State VM;
    if(argc<=1)
    {
        throw exception("need a config.lua");
    }

    lua_State *L = VM();
    Lua::Config::DoFile(L,argv[1]);

    Profile pro(L);

    {
        ios::ocstream fp("profile.dat",false);
        const size_t NP = 100;
        for(size_t i=0;i<=NP;++i)
        {
            const double z = double(i)/NP;
            fp("%g %g %g\n", z, pro.width(z), pro.height(z));
        }
    }

    {
        ios::ocstream fp("rmax.dat",false);
        for(size_t i=1;i<=pro.zarr.size();++i)
        {
            fp("%g %g\n", pro.zarr[i], pro.rmax[i]);
        }
    }
}
YOCTO_PROGRAM_END()
