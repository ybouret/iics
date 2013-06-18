#include "simulation.hpp"

Simulation:: ~Simulation() throw()
{
    
    
}


Simulation:: Simulation(const mpi   &MPI,
                        const Coord  N,
                        const Vertex Q) :
Workspace(MPI,N,Q),
VisIt::Simulation(MPI),
ftol(1e-4)
{
    
}

void Simulation:: initialize()
{
    validate_bubbles(MPI);
    if(!is_valid)
    {
        done = true;
        return;
    }
    
    broadcast_bubbles(MPI);
    segment();
    
    compute_pressure(MPI,ftol);
}



////////////////////////////////////////////////////////////////////////////////
//
// Specific perform
//
////////////////////////////////////////////////////////////////////////////////
#include "yocto/string/conv.hpp"
#include "../shape.hpp"


void Simulation:: perform( const string &cmd, const array<string> &args)
{
    if( cmd == "raz" )
    {
        const char *kind = 0;
        if( args.size() > 0 )
            kind = args[1].c_str();
        init_one_bubble(kind);
        initialize();
        return;
    }
    
    if(cmd == "gamma" )
    {
        
        bubbles.gamma = args.size() >= 1 ? strconv::to<Real>(args[1],"gamma") : 0;
        MPI.Printf0(stderr, "Changing gamma to %g\n", bubbles.gamma);
        initialize();
        return;
    }

    
    if(cmd == "rot" )
    {
        if( args.size() >= 1 )
        {
            
            if(MPI.IsFirst)
            {
                const Real angle = strconv::to<Real>(args[1],"angle");
                for(Bubble *b=bubbles.head;b;b=b->next)
                {
                    Shape::Rotate(b, angle*numeric<Real>::pi/180.0);
                }
            }
            initialize();
        }
        return;
    }

    if(cmd == "grow" )
    {
        if( args.size() >= 1 )
        {
            
            if(MPI.IsFirst)
            {
                const Real factor = strconv::to<Real>(args[1],"factor");
                for(Bubble *b=bubbles.head;b;b=b->next)
                {
                    Shape::Grow(b, factor);
                }
            }
            initialize();
        }
        return;
    }

    
#if 0
    if(cmd == "dx")
    {
        if( args.size() >= 1 )
        {
            if(MPI.IsFirst)
            {
                const Real dx = strconv::to<Real>(args[1],"dx");
                const Vertex v(dx,0);
                for(Bubble *b=bubbles.head;b;b=b->next)
                {
                    Shape::Move(b,v);
                }
            }
            fast_update();
        }
    }

    
    
    
    if(cmd == "dy")
    {
        if( args.size() >= 1 )
        {
            if(MPI.IsFirst)
            {
                const Real dy = strconv::to<Real>(args[1],"dy");
                const Vertex v(0,dy);
                for(Bubble *b=bubbles.head;b;b=b->next)
                {
                    Shape::Move(b,v);
                }
            }
            fast_update();
        }
    }
    
    
    
       
    
    if(cmd == "rot" )
    {
        if( args.size() >= 1 )
        {
            
            if(MPI.IsFirst)
            {
                const Real angle = strconv::to<Real>(args[1],"angle");
                for(Bubble *b=bubbles.head;b;b=b->next)
                {
                    Shape::Rotate(b, angle*numeric<Real>::pi/180.0);
                }
            }
            fast_update();
        }
    }
    

    
    if(cmd == "rb" )
    {
        size_t n = args.size() >= 1 ? strconv::to<size_t>(args[1],"niter") : 1;
        DeltaP.ldz();
        while(n-->0)
        {
            const int cvg = update_pressure(MPI, Red, ftol) & update_pressure(MPI, Black, ftol);
            MPI.Printf(stderr, "Converged= %d\n", cvg);
        }
    }
    
    if(cmd == "gamma" )
    {
        
        bubbles.gamma = args.size() >= 1 ? strconv::to<Real>(args[1],"gamma") : 0;
        segment();
        compute_pressure(MPI, ftol);
        MPI.Printf(stderr, "Pressure OK\n");
    }
    
    if(cmd == "solve" )
    {
        compute_pressure(MPI, ftol);
        MPI.Printf(stderr, "Pressure OK\n");
    }
    
    if(cmd == "save" )
    {
        if(MPI.IsFirst )
        {
            for( const Bubble *b = bubbles.head;b;b=b->next)
            {
                const string pfx = vformat("b%u", unsigned(b->UID) );
                b->save_all(pfx);
            }
            junctions.save_all( "j" + MPI.CommWorldID);
            save_markers();
        }
    }
    
       
    if( cmd == "zp" )
    {
        P.ldz();
        pressurize_bubbles();
        pressurize_contours();
        compute_gradP(MPI);
        compute_laplacian();
    }
    
    if( cmd == "full" )
    {
        if(MPI.IsFirst)
        {
            for(Bubble *b=bubbles.head;b;b=b->next)
            {
                Shape::Rotate(b, 10*numeric<Real>::pi/180.0);
            }
        }
        validate_bubbles(MPI);
        if( !is_valid)
            throw exception("Invalid bubbles");
        broadcast_bubbles(MPI);
        segment();
        compute_pressure(MPI, 1e-4);
        
    }
#endif
    
}
