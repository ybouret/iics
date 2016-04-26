#include "yocto/graphics/image/png.hpp"
#include "yocto/graphics/image/jpeg.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/graphics/rawpix.hpp"

#include "yocto/program.hpp"

using namespace yocto;
using namespace graphics;

YOCTO_PROGRAM_START()
{

    //__________________________________________________________________________
    //
    // will use some image format
    //__________________________________________________________________________
    image &IMG = image::instance();
    IMG.declare( new png_format()   );
    IMG.declare( new jpeg_format()  );

    //__________________________________________________________________________
    //
    // will use Lua
    //__________________________________________________________________________
    Lua::State VM;
    lua_State *L = VM();

    if(argc>1)
    {
        Lua::Config::DoFile(L, argv[1]);
        for(int i=2;i<argc;++i)
        {
            Lua::Config::DoString(L,argv[i]);
        }
    }

    const unit_t w = unit_t(Lua::Config::Get<lua_Number>(L, "w"));
    const unit_t h = unit_t(Lua::Config::Get<lua_Number>(L, "h"));

    pixmap1 img1(w,h);
    pixmap1 img2(w,h);


}
YOCTO_PROGRAM_END()
