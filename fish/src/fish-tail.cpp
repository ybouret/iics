#include "fish.hpp"
#include "yocto/exception.hpp"

void Fish:: generateTail( double Zmax, size_t N )
{

    // clean up
    clear();

    if(Zmax<=0||Zmax>0.5)
        throw exception("Invalid Zmax=%g",Zmax);

    zfunction<double> rfn( this, & Profile::getZ, Zmax);
    zfind<double>     solve(0);
    const double      rho_max =  solve(rfn.call,0,1);
    std::cerr << "rho_max=" << rho_max << std::endl;
    N = max_of<size_t>(1,N);

    const double delta = (1-rho_max) / N;
    const size_t n = max_of<size_t>(2,ceil(maxP/delta))-2;
    const size_t M = 2+2*n;
    std::cerr << "\tn=" << n << ", M=" << M << std::endl;


    // compute the slices position
    for(size_t i=1;i<=N;++i)
    {
        const double ratio = rho_max + (i-1)*delta;
        const pSlice pS( new Slice( getZ(ratio) ) );
        slices.push_back(pS);
    }



    // head point
    pPoint p0( new Point() );
    points.push_back(p0);

    p0->r.z = slices.front()->z;

    // body points
    for(size_t i=1;i<=N;++i)
    {

        Slice &slice = *slices[i];

        const double w = W(slice.z);
        const double h = H(slice.z);

        //std::cerr << "w=" << width << ", h=" << height << std::endl;
        for(size_t j=0;j<M;++j)
        {
            const double theta = (j*numeric<double>::two_pi)/M;
            pPoint pp( new Point() );

            pp->r.x = w  * cos(theta);
            pp->r.y = h * sin(theta);
            pp->r.z = slice.z;

            slice.points.push_back(pp);
            points.push_back(pp);
        }

    }


    // tail points
    pPoint pN( new Point() );
    points.push_back(pN);
    pN->r.z = 1;

    std::cerr << "#points=" << points.size() << std::endl;

    // generate triangles

    std::cerr << "-- Computing Triangles" << std::endl;
    std::cerr << "\t Head..." << std::endl;

    // head
    {
        const Slice &slice = *slices[1];
        for(size_t i=1;i<=M;++i)
        {
            size_t   ip = i+1;
            if(ip>M) ip = 1;
            Triangle tr(p0,slice.points[i],slice.points[ip]);
            if(tr.n.z>=0) tr.inverse();
            triangles.push_back(tr);
        }
    }

    std::cerr << "\t Body..." << std::endl;
    // inside
    for(size_t j=1;j<N;++j)
    {
        const array<pPoint> &P0 = slices[j]->points;
        const array<pPoint> &P1 = slices[j+1]->points;

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
            }

            {
                const Triangle tr(P00,P10,P11);
                triangles.push_back(tr);
            }


        }
    }

    std::cerr << "\t Tail..." << std::endl;
    // tail
    {
        const Slice &slice = *slices[N];
        for(size_t i=1;i<=M;++i)
        {
            size_t   ip = i+1;
            if(ip>M) ip = 1;
            const Triangle tr(pN,slice.points[i],slice.points[ip]);
            triangles.push_back(tr);
        }
    }
    
    std::cerr << "#triangles=" << triangles.size() << std::endl;
    
    
}
