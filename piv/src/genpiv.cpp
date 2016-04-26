#include "yocto/graphics/image/png.hpp"
#include "yocto/graphics/image/jpeg.hpp"
#include "yocto/lua/lua-state.hpp"
#include "yocto/lua/lua-config.hpp"
#include "yocto/graphics/rawpix.hpp"
#include "yocto/random/bits.hpp"
#include "yocto/random/default.hpp"
#include "yocto/code/rand.hpp"
#include "yocto/random/gaussian.hpp"
#include "yocto/code/utils.hpp"
#include "yocto/program.hpp"

using namespace yocto;
using namespace graphics;

YOCTO_PROGRAM_START()
{

    alea_init();

    //__________________________________________________________________________
    //
    // will use some image format
    //__________________________________________________________________________
    image &IMGDB = image::instance();
    IMGDB.declare( new png_format()   );
    IMGDB.declare( new jpeg_format()  );
    imageIO &IMG = IMGDB;

    //__________________________________________________________________________
    //
    // will use Lua: one file and some overrides
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

    //__________________________________________________________________________
    //
    // reading parameters
    //__________________________________________________________________________
    const unit_t w = unit_t(Lua::Config::Get<lua_Number>(L, "w"));
    const unit_t h = unit_t(Lua::Config::Get<lua_Number>(L, "h"));
    const double intensity_rho     = double(Lua::Config::Get<lua_Number>(L,"intensity_rho"));
    const float  intensity_mu      = float(Lua::Config::Get<lua_Number>(L,"intensity_mu"));
    const float  intensity_sigma   = float(Lua::Config::Get<lua_Number>(L,"intensity_sigma"));
    const size_t particle_max_size = max_of<size_t>(1,size_t(Lua::Config::Get<lua_Number>(L,"particle_max_size")));

    //__________________________________________________________________________
    //
    // target images
    //__________________________________________________________________________
    pixmap1 img1(w,h);
    pixmap1 img2(w,h);

    //__________________________________________________________________________
    //
    // source image
    //__________________________________________________________________________
    Random::Default  ran( Random::SimpleBits() );
    Random::Gaussian gran(ran);
    
    const unit_t W = w+w;
    const unit_t H = h+h;
    pixmap1      src(W,H);
    for(unit_t j=0;j<H;++j)
    {
        for(unit_t i=0;i<W;++i)
        {
            const double r = ran();
            //std::cerr << "r=" << r << std::endl;
            if(r>intensity_rho)
            {
                continue;
            }

            const unit_t  s = unit_t(1+alea_lt(particle_max_size));
            for(unit_t jj=0;jj<s;++jj)
            {
                const unit_t J = j+jj;
                if(J<0||J>=H) continue;
                for(unit_t ii=0;ii<s;++ii)
                {
                    const unit_t  I = i + ii;
                    if(I<0||I>=W)
                    {
                        continue;
                    }

                    const float   q = clamp<float>(0,intensity_mu+intensity_sigma*gran(),1);
                    const uint8_t u = gist::float2byte(q);
                    src[J][I] = u;
                }
            }
        }
    }

    //__________________________________________________________________________
    //
    // img1
    //__________________________________________________________________________
    unit_t x=w/2;
    unit_t y=h/2;
    for(unit_t j=0;j<h;++j)
    {
        for(unit_t i=0;i<w;++i)
        {
            img1[j][i]=src[j+y][i+x];
        }
    }

    


    //__________________________________________________________________________
    //
    // finally: save
    //__________________________________________________________________________
    //IMG.save("src.png",  src,  NULL);
    IMG.save("img1.png", img1, NULL);

    for(int k=2;k<=20;++k)
    {
        img2.ldz();
        ++x;
        for(unit_t j=0;j<h;++j)
        {
            for(unit_t i=0;i<w;++i)
            {
                const float f = gist::unit_float[src[j+y][i+x]];
                img2[j][i] = gist::float2byte(clamp<float>(0,f+0.5f*(alea<float>()-0.5f),1));
            }
        }
        IMG.save( vformat("img%d.png",k), img2, NULL);

    }


}
YOCTO_PROGRAM_END()
