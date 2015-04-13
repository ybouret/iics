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

    Fish fish(L);

    {
        ios::ocstream fp("profile.dat",false);
        const size_t NP = 100;
        for(size_t i=0;i<=NP;++i)
        {
            const double z = double(i)/NP;
            fp("%g %g %g\n", z, fish.W(z), fish.H(z));
        }
    }

    {
        ios::ocstream fp("rmax.dat",false);
        for(size_t i=1;i<=fish.zarr.size();++i)
        {
            fp("%g %g %g %g\n", fish.zarr[i], fish.rmax[i], fish.arcL[i], fish.computePerimeter(fish.zarr[i]));
        }
    }

    {
        ios::ocstream fp("zpos.dat",false);
        for(size_t i=0;i<=100;++i)
        {
            const double ratio = i/100.0;
            fp("%g %g\n", ratio, fish.getZ(ratio));
        }
    }



    fish.generateShell( Lua::Config::Get<lua_Number>(L,"N") );
    fish.centerAndRescaleBy( Lua::Config::Get<lua_Number>(L,"Length"));
    fish.save_vtk("fish_shell.vtk");


    fish.generateHead(0.33, 10, 0.05);
    fish.save_vtk("fish_head.vtk");

}
YOCTO_PROGRAM_END()
