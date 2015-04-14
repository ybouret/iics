#include "fish.hpp"
#include "yocto/exception.hpp"

static inline
double ReverseMax(const double               U,
                  numeric<double>::function &lambda,
                  double                     z,
                  double                     step )
{
    assert(z>0);
    assert(z<1);
    assert(U<lambda(z));
    step = Fabs(step);
    const double z_ini = z;

    for(z-=step;z>0;z-=step)
    {
        if(lambda(z)<=U)
            break;
    }

    return max_of<double>(z_ini + (z-z_ini)/2,0.0);
}

static inline
double ForwardMax(const double               U,
                  numeric<double>::function &lambda,
                  double                     z,
                  double                     step )
{
    assert(z>0);
    assert(z<1);
    assert(U<lambda(z));
    step = Fabs(step);
    const double z_ini = z;

    for(z+=step;z<1;z+=step)
    {
        //std::cerr << "fwd z=" << z << std::endl;
        if(lambda(z)<=U)
            break;
    }

    return min_of<double>(z_ini + (z-z_ini)/2,1.0);
}


void Fish:: generateJunction( double Zmax, double thickness)
{
    if(Zmax<=0||Zmax>0.5)
        throw exception("Invalid Zmax=%g",Zmax);


    const double w = W(Zmax);
    const double h = W(Zmax);

    const double alpha = 1.0 - thickness/max_of(w,h);

    const double w_in  = w * alpha;
    const double h_in  = h * alpha;

    const double w_out = w_in + thickness/2;
    const double h_out = h_in + thickness/2;

    const double zwr   = ReverseMax(w_out,W,Zmax,1e-5);
    const double hwr   = ReverseMax(h_out,H,Zmax,1e-5);

    std::cerr << "w_in =" << w_in  << ", h_in=" << h_in << std::endl;
    std::cerr << "w_out=" << w_out << ", h_out=" << h_out << std::endl;

    //std::cerr << "max_reverse_zW=" << zwr << std::endl; std::cerr << "max_reverse_hW=" << hwr << std::endl;

    const double reverse_z = max_of(zwr,hwr);
    std::cerr << "reverse_z=" << reverse_z << std::endl;

    const double zwf = ForwardMax(w_out, W, Zmax, 1e-5);
    const double hwf = ForwardMax(h_out, H, Zmax, 1e-5);

    //std::cerr << "max_forward_zW=" << zwf << std::endl; std::cerr << "max_forward_hW=" << hwf << std::endl;

    const double forward_z = min_of(zwf,hwf);
    std::cerr << "forward_z=" << forward_z << std::endl;

}