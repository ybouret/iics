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


    const size_t N         = Lua::Config::Get<lua_Number>(L,"N");
    const double Length    = Lua::Config::Get<lua_Number>(L,"Length");
    const double thickness = Lua::Config::Get<lua_Number>(L,"thickness");
    const double head_size = clamp<double>(0.1,Lua::Config::Get<lua_Number>(L,"head_size"),0.9);
    const double junc_size = Lua::Config::Get<lua_Number>(L,"junc_size");
    
    fish.generateShell(  N );
    fish.centerAndRescale(head_size,Length);

    fish.save_vtk("fish_shell.vtk");


   

    const size_t NH = ceil(head_size*N);
    fish.generateHead(head_size, NH, thickness);
    fish.centerAndRescale(head_size,Length);
    fish.save_vtk("fish_head.vtk");
    fish.save_stl("fish_head.stl");

    const size_t NT = ceil((1.0-head_size)*N);

    fish.generateTail(head_size, NT );
    fish.centerAndRescale(head_size,Length);
    fish.save_vtk("fish_tail.vtk");
    fish.save_stl("fish_tail.stl");

    fish.generateJunction(head_size,thickness,junc_size);
    fish.centerAndRescale(head_size,Length);
    fish.save_vtk("fish_junc.vtk");
    fish.save_stl("fish_junc.stl");


}
YOCTO_PROGRAM_END()
