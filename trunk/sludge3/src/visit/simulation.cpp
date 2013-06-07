#include "simulation.hpp"

Simulation:: ~Simulation() throw()
{
    
    
}


Simulation:: Simulation(const mpi   &MPI,
                        const Coord  N,
                        const Vertex Q) :
Workspace(MPI,N,Q),
VisIt::Simulation(MPI)
{
    
}

////////////////////////////////////////////////////////////////////////////////
//
// Get Meta Data
//
////////////////////////////////////////////////////////////////////////////////

void Simulation:: get_meta_data(visit_handle &md) const
{
    
    //! append an UI generic command
    add_generic_command("raz", md);
    
    //! append the grid
    VisIt_SimulationMetaData_addMesh(md, mesh_meta_data(mesh, "grid", par_size));
    
    
    //! append P on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Real>("P", "grid"));
    
    
    //! append B on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Real>("B", "grid"));
    
    //! append gradP on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("gradP", "grid"));
    
    //! append Enter on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("E1", "grid"));
    
    //! append Leave on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("L1", "grid"));
    
    //! append Enter on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("E2", "grid"));
    
    //! append Leave on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Vertex>("L2", "grid"));
    
    
    //! append DeltaP on grid
    VisIt_SimulationMetaData_addVariable(md, variable_meta_data<Real>("DeltaP", "grid"));
    
    //! append bubbles
    for( const Bubble *b = bubbles.head;b;b=b->next)
    {
        const string bn = vformat("bubble%u", unsigned(b->UID));
        visit_handle h = VISIT_INVALID_HANDLE;
        if(VisIt_CurveMetaData_alloc(&h) == VISIT_OKAY )
        {
            VisIt_CurveMetaData_setName(h, bn.c_str());
            //VisIt_CurveMetaData_setXLabel(cmd, "Angle");
            ////VisIt_CurveMetaData_setXUnits(cmd, "radians");
            //VisIt_CurveMetaData_setYLabel(cmd, "Amplitude");
            VisIt_SimulationMetaData_addCurve(md, h);
        }
    }
    
    //! append junctions
    if( !parallel )
    {
        visit_handle h = VISIT_INVALID_HANDLE;
        if(VisIt_CurveMetaData_alloc(&h) == VISIT_OKAY )
        {
            VisIt_CurveMetaData_setName(h, "junctions" );
            VisIt_SimulationMetaData_addCurve(md, h);
        }
    }
    
}

////////////////////////////////////////////////////////////////////////////////
//
// Get Mesh
//
////////////////////////////////////////////////////////////////////////////////

visit_handle Simulation:: get_mesh( int domain, const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if( name == "grid" )
    {
        if( VisIt_RectilinearMesh_alloc(&h) == VISIT_OKAY )
        {
            assert( h != VISIT_INVALID_HANDLE );
            
            visit_handle   hx = VISIT_INVALID_HANDLE;
            visit_handle   hy = VISIT_INVALID_HANDLE;
            VisIt_VariableData_alloc(&hx);
            VisIt_VariableData_setDataD(hx, VISIT_OWNER_SIM, 1, X.items, X.entry );
            VisIt_VariableData_alloc(&hy);
            VisIt_VariableData_setDataD(hy, VISIT_OWNER_SIM, 1, Y.items, Y.entry );
            
            VisIt_RectilinearMesh_setCoordsXY(h, hx, hy);
            
#if 0
            int minRealIndex[3] = { 1,         2 ,        0 };
            int maxRealIndex[3] = { int(X.width)-2, int(Y.width)-3, 0 };
            
            if(  par_rank < par_size - 1)
                maxRealIndex[1]++;
            
            VisIt_RectilinearMesh_setRealIndices(h, minRealIndex, maxRealIndex);
#endif
        }
    }
    return h;
    
}


////////////////////////////////////////////////////////////////////////////////
//
// Get Variables
//
////////////////////////////////////////////////////////////////////////////////

visit_handle Simulation:: get_variable( int domain, const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    
    if( name == "P" )
    {
        const int nComponents= 1;
        const int nTuples    = P.items;
        assert(P.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, P.entry);
        }
    }
    
    if( name == "B" )
    {
        const int nComponents= 1;
        const int nTuples    = B.items;
        assert(B.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, B.entry);
        }
    }
    
    if( name == "gradP" )
    {
        const int nComponents= 2;
        const int nTuples    = gradP.items;
        assert(gradP.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real *)(gradP.entry));
        }
    }
    
    if( name == "E1" )
    {
        const int nComponents= 2;
        const int nTuples    = E1.items;
        assert(E1.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real *)(E1.entry));
        }
    }
    
    if( name == "L1" )
    {
        const int nComponents= 2;
        const int nTuples    = L1.items;
        assert(L1.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real *)(L1.entry));
        }
    }
    
    if( name == "E2" )
    {
        const int nComponents= 2;
        const int nTuples    = E1.items;
        assert(E1.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real *)(E2.entry));
        }
    }
    
    if( name == "L2" )
    {
        const int nComponents= 2;
        const int nTuples    = L1.items;
        assert(L1.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (Real *)(L2.entry));
        }
    }
    
    
    if( name == "DeltaP" )
    {
        const int nComponents= 1;
        const int nTuples    = DeltaP.items;
        assert(DeltaP.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, DeltaP.entry);
        }
    }
    
    
    return h;
    
}

visit_handle Simulation:: get_curve( const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if( !parallel )
    {
        if( name == "junctions")
        {
            if(VisIt_CurveData_alloc(&h) != VISIT_ERROR)
            {
                visit_handle hxc, hyc;
                VisIt_VariableData_alloc(&hxc);
                VisIt_VariableData_alloc(&hyc);
                
                const size_t nj = junctions.count_all();
                vector<Real> jx(nj,0);
                vector<Real> jy(nj,0);
                junctions.to_curve(jx, jy);
                
                VisIt_VariableData_setDataD(hxc,VISIT_OWNER_COPY,1,nj,&jx[1]);
                VisIt_VariableData_setDataD(hyc,VISIT_OWNER_COPY,1,nj,&jy[1]);
                VisIt_CurveData_setCoordsXY(h, hxc, hyc);
                
                return h;
            }
        }
    }
    
    for( const Bubble *b = bubbles.head;b;b=b->next)
    {
        const string bn = vformat("bubble%u", unsigned(b->UID));
        if( bn == name )
        {
            if(VisIt_CurveData_alloc(&h) != VISIT_ERROR)
            {
                visit_handle hxc, hyc;
                VisIt_VariableData_alloc(&hxc);
                VisIt_VariableData_alloc(&hyc);
                
                const size_t np = b->size+1;
                vector<Real> bx(np,0);
                vector<Real> by(np,0);
                const Tracer *tr = b->root;
                for(size_t i=1;i<=np;++i,tr=tr->next)
                {
                    bx[i] = tr->pos.x;
                    by[i] = tr->pos.y;
                }
                VisIt_VariableData_setDataD(hxc,VISIT_OWNER_COPY,1,np,&bx[1]);
                VisIt_VariableData_setDataD(hyc,VISIT_OWNER_COPY,1,np,&by[1]);
                VisIt_CurveData_setCoordsXY(h, hxc, hyc);
            }
        }
    }
    return h;
    
}



void Simulation:: fast_update()
{
    validate_bubbles(MPI);
    if(!is_valid)
    {
        done = true;
        return;
    }
    broadcast_bubbles(MPI);
    segment();
    P.ldz();
    pressurize_bubbles();
    pressurize_contours();
    compute_gradP(MPI);
    
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
    }
    
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
            fast_update();
        }
    }
    
    
    const Real ftol = 1e-5;
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
        }
    }
}
