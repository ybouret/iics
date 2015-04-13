#include "fish.hpp"
#include "yocto/exception.hpp"


void Fish:: generateHead(  double Zmax, size_t N, double thickness)
{
    //--------------------------------------------------------------------------
    //
    // Setup
    //
    //--------------------------------------------------------------------------

    clear();
    if(Zmax<=0||Zmax>0.5)
        throw exception("Invalid Zmax=%g",Zmax);

    zfunction<double> rfn( this, & Profile::getZ, Zmax);
    zfind<double>     solve(0);
    const double      rho_max =  solve(rfn.call,0,1);
    std::cerr << "rho_max=" << rho_max << std::endl;
    N = max_of<size_t>(1,N);

    const double delta = rho_max / (N);
    const size_t n = max_of<size_t>(2,ceil(maxP/delta))-2;
    const size_t M = 2+2*n;
    std::cerr << "\tn=" << n << ", M=" << M << std::endl;

    //--------------------------------------------------------------------------
    //
    // Outer Shell
    //
    //--------------------------------------------------------------------------
    {
        //______________________________________________________________________
        //
        // compute the slices position
        //______________________________________________________________________
        for(size_t i=1;i<=N;++i)
        {
            const double ratio = i*delta;
            const pSlice pS( new Slice( getZ(ratio) ) );
            slices.push_back(pS);
        }

        //______________________________________________________________________
        //
        // compute the number of points per slice M = 2*n+2;
        //______________________________________________________________________



        //______________________________________________________________________
        //
        // head point
        //______________________________________________________________________
        pPoint p0( new Point() );
        points.push_back(p0);

        //______________________________________________________________________
        //
        // body points
        //______________________________________________________________________
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

                pp->r.x = w * cos(theta);
                pp->r.y = h * sin(theta);
                pp->r.z = slice.z;

                slice.points.push_back(pp);
                points.push_back(pp);
            }

        }

        std::cerr << "-- Computing Triangles" << std::endl;
        std::cerr << "\t Head..." << std::endl;

        //______________________________________________________________________
        //
        // head
        //______________________________________________________________________
        {
            const Slice &slice = *slices[1];
            for(size_t i=1;i<=M;++i)
            {
                size_t   ip = i+1;
                if(ip>M) ip = 1;
                const Triangle tr(p0,slice.points[i],slice.points[ip]);
                triangles.push_back(tr);
            }
        }

        std::cerr << "\t Body..." << std::endl;
        //______________________________________________________________________
        //
        // inside
        //______________________________________________________________________
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

    }

    //--------------------------------------------------------------------------
    //
    // Inner Shell
    //
    //--------------------------------------------------------------------------
    {
        vector<pPoint>   inner_points;
        vector<pSlice>   inner_slices;
        vector<Triangle> inner_tr;

        pPoint p0( new Point() );
        inner_points.push_back(p0);

        rfn.target = thickness;
        const double rho_min = solve(rfn.call,0,1);
        p0->r.z = thickness;

        const double delta_in = (rho_max - rho_min)/N;
        //______________________________________________________________________
        //
        // inner slices
        //______________________________________________________________________
        for(size_t i=1;i<=N;++i)
        {
            const double ratio = rho_min + (i*delta_in);
            const pSlice pS( new Slice( getZ(ratio) ) );
            inner_slices.push_back(pS);
        }


        //______________________________________________________________________
        //
        // inner body points
        //______________________________________________________________________
        for(size_t i=1;i<=N;++i)
        {

            Slice &slice = *inner_slices[i];

            double w = W(slice.z);
            double h = H(slice.z);
            const double r = Hypotenuse(w, h);
            const double fac = (r-thickness)/r;

            w *= fac;
            h *= fac;

            //std::cerr << "w=" << width << ", h=" << height << std::endl;
            for(size_t j=0;j<M;++j)
            {
                const double theta = (j*numeric<double>::two_pi)/M;
                pPoint pp( new Point() );

                pp->r.x = w * cos(theta);
                pp->r.y = h * sin(theta);
                pp->r.z = slice.z;

                slice.points.push_back(pp);
                inner_points.push_back(pp);
            }

        }

        std::cerr << "-- Computing Inner Triangles" << std::endl;
        std::cerr << "\t Head..." << std::endl;
        //______________________________________________________________________
        //
        // inner triangles head
        //______________________________________________________________________
        {
            const Slice &slice = *inner_slices[1];
            for(size_t i=1;i<=M;++i)
            {
                size_t   ip = i+1;
                if(ip>M) ip = 1;
                const Triangle tr(p0,slice.points[i],slice.points[ip]);
                inner_tr.push_back(tr);
            }
        }

        std::cerr << "\t Body..." << std::endl;
        //______________________________________________________________________
        //
        // inside
        //______________________________________________________________________
        for(size_t j=1;j<N;++j)
        {
            const array<pPoint> &P0 = inner_slices[j]->points;
            const array<pPoint> &P1 = inner_slices[j+1]->points;

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
                    inner_tr.push_back(tr);
                }

                {
                    const Triangle tr(P00,P10,P11);
                    inner_tr.push_back(tr);
                }
            }
        }

        //______________________________________________________________________
        //
        // fusion
        //______________________________________________________________________
        for(size_t i=1;i<=inner_points.size();++i)
        {
            points.push_back(inner_points[i]);
        }

        for(size_t i=1;i<=inner_tr.size();++i)
        {
            inner_tr[i].inverse();
            triangles.push_back( inner_tr[i] );
        }

        //______________________________________________________________________
        //
        // closing the shell
        //______________________________________________________________________
        std::cerr << "-- Closing the shell" << std::endl;
        const Slice &S_out = *slices.back();
        const Slice &S_ins  = *inner_slices.back();
        {
            const array<pPoint> &P0 = S_out.points;
            const array<pPoint> &P1 = S_ins.points;

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
        
    }
    
}
