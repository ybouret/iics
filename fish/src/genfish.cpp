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
            fp("%g %g %g\n", z, pro.W(z), pro.H(z));
        }
    }

    {
        ios::ocstream fp("rmax.dat",false);
        for(size_t i=1;i<=pro.zarr.size();++i)
        {
            fp("%g %g %g %g\n", pro.zarr[i], pro.rmax[i], pro.arcL[i], pro.computePerimeter(pro.zarr[i]));
        }
    }

    {
        ios::ocstream fp("zpos.dat",false);
        for(size_t i=0;i<=100;++i)
        {
            const double ratio = i/100.0;
            fp("%g %g\n", ratio, pro.getZ(ratio));
        }
    }

}
YOCTO_PROGRAM_END()
