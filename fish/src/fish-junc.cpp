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

static inline void GenTr(vector<Triangle>    &triangles,
                         const array<pPoint> &P0,
                         const array<pPoint> &P1,
                         const bool           do_inv)
{
    assert(P0.size()==P1.size());
    const size_t M = P0.size();

    // loop over quads
    for(size_t i=1;i<=M;++i)
    {
        size_t   ip = i+1;
        if(ip>M) ip = 1;
        const pPoint  &P00 = P0[i];
        const pPoint  &P01 = P0[ip];
        const pPoint  &P10 = P1[i];
        const pPoint  &P11 = P1[ip];

        {
            const Triangle tr(P00,P01,P11);
            triangles.push_back(tr);
            if(do_inv) triangles.back().inverse();
        }

        {
            const Triangle tr(P00,P10,P11);
            triangles.push_back(tr);
            if(do_inv) triangles.back().inverse();
        }


    }

}



void Fish:: generateJunction( double Zmax, double thickness, double junc_size )
{
    clear();

    if(Zmax<=0||Zmax>0.5)
        throw exception("Invalid Zmax=%g",Zmax);



    const double w = W(Zmax);
    const double h = H(Zmax);

    std::cerr << "w,h= " << w << ", " << h << std::endl;

    const double w_in  = w - thickness/2;
    const double h_in  = h - thickness/2;
    std::cerr << "w_in,h_in= " << w_in << ", " << h_in << std::endl;


    const double w_out = w_in + thickness/4;
    const double h_out = h_in + thickness/4;

    std::cerr << "w_out,h_out= " << w_out << ", " << h_out << std::endl;

    const double zwr   = ReverseMax(w_out,W,Zmax,1e-5);
    const double hwr   = ReverseMax(h_out,H,Zmax,1e-5);

    std::cerr << "w_in =" << w_in  << ", h_in=" << h_in << std::endl;
    std::cerr << "w_out=" << w_out << ", h_out=" << h_out << std::endl;

    //std::cerr << "max_reverse_zW=" << zwr << std::endl; std::cerr << "max_reverse_hW=" << hwr << std::endl;

    const double reverse_z = max_of(max_of(zwr,hwr),Zmax-Fabs(junc_size)/2);

    std::cerr << "reverse_z=" << reverse_z << std::endl;

    const double zwf = ForwardMax(w_out, W, Zmax, 1e-5);
    const double hwf = ForwardMax(h_out, H, Zmax, 1e-5);

    //std::cerr << "max_forward_zW=" << zwf << std::endl; std::cerr << "max_forward_hW=" << hwf << std::endl;

    const double forward_z = min_of(min_of(zwf,hwf),Zmax+Fabs(junc_size)/2);
    std::cerr << "forward_z=" << forward_z << std::endl;

    Slice RO(reverse_z);
    Slice RI(reverse_z);

    Slice FO(reverse_z);
    Slice FI(reverse_z);


    // build slices
    size_t M = 60;
    for(size_t j=0;j<M;++j)
    {
        const double theta = (j*numeric<double>::two_pi)/M;
        const double ct    = cos(theta);
        const double st    = sin(theta);

        pPoint pRO( new Point() );  RO.points.push_back(pRO);
        pPoint pRI( new Point() );  RI.points.push_back(pRI);

        pPoint pFO( new Point() );  FO.points.push_back(pFO);
        pPoint pFI( new Point() );  FI.points.push_back(pFI);

        pRO->r.z = reverse_z;
        pRI->r.z = reverse_z;

        pFO->r.z = forward_z;
        pFI->r.z = forward_z;

        pRO->r.x = pFO->r.x = w_out * ct;
        pRO->r.y = pFO->r.y = h_out * st;

        pRI->r.x = pFI->r.x = w_in * ct;
        pRI->r.y = pFI->r.y = h_in * st;



        points.push_back(pRO);
        points.push_back(pRI);
        points.push_back(pFO);
        points.push_back(pFI);

    }

    // build triangles
    GenTr(triangles, RO.points, FO.points, false);
    GenTr(triangles, RI.points, FI.points, true);
    
    // close Reverse Side
    {
        vector<Triangle> tr;
        GenTr(tr, RO.points, RI.points, false);
        for(size_t i=1;i<=tr.size();++i)
        {
            if(tr[i].n.z>=0) tr[i].inverse();
            triangles.push_back(tr[i]);
        }
    }

    // close Forward Side
    {
        vector<Triangle> tr;
        GenTr(tr, FO.points, FI.points, false);
        for(size_t i=1;i<=tr.size();++i)
        {
            if(tr[i].n.z<=0) tr[i].inverse();
            triangles.push_back(tr[i]);
        }
    }



}