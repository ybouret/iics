#include "yocto/utest/run.hpp"
#include "yocto/code/rand.hpp"
#include "shape.hpp"


YOCTO_UNIT_TEST_IMPL(ctrl)
{
    Real   lam = 0.2;
    Real   mu  = 0.15;
    Real   gam = 1;
    Bubble bubble(lam,gam,1);
    
    std::cerr << "Circle" << std::endl;
    Shape::Circle(&bubble, Vertex(0,0), 4);
    bubble.save_dat( "c0.dat");
    
    bubble.reduce(mu);
    bubble.save_dat( "c1.dat");
    
    
    std::cerr << "Ellipse" << std::endl;
    Shape::Ellipse(&bubble, Vertex(0,0), Vertex(5,2));
    bubble.save_dat( "e0.dat");
    
    bubble.reduce(mu);
    bubble.save_dat( "e1.dat");
    
    std::cerr << "Blob" << std::endl;
    Shape::Blob(&bubble, Vertex(0,0), 4, 0.5, 0.0);
    bubble.save_dat( "b0.dat");
    
    for(size_t i=0;i<10;++i)
        bubble.reduce1(mu);
    bubble.save_dat( "b1.dat");
}
YOCTO_UNIT_TEST_DONE()
