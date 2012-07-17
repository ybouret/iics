#include "simulation.hpp"
#include "yocto/swamp/mpi.hpp"

Simulation:: Simulation(unit_t     Nx, 
                        unit_t     Ny,
                        Real       Lx,
                        Real       Ly,
                        const mpi &ref) :
Cell(Nx,Ny,Lx,Ly,ref),
VisIt::Simulation(ref),
requests( num_requests())
{
    prepare_ghosts();
}

Simulation:: ~Simulation() throw()
{
    
}


void Simulation:: initialize()
{
    U.ldz();
    P.ldz();
    bubbles.none();
    for( unit_t j=lower.y;j<=upper.y;++j)
    {
        for( unit_t i=lower.x;i<=upper.x;++i)
        {
            P[j][i]   = Y[j] / Length.y;
            U[j][i].y = 0.04;
        }
    }
    if( master )
    {
        Bubble *b = bubbles.create();
        b->lambda = Length.x / (width.x*2);
        b->map_circle( V2D(Length.x/2,0), 0.2 * Length.y);
    }
}

void Simulation:: init_exchange()
{
    _mpi::init_exchange(MPI, *this, requests);
}

void Simulation:: wait_exchange()
{
    _mpi::wait_exchange(MPI, *this, requests);
}

void Simulation:: check_and_dispatch_bubbles()
{
    if( master )
    {
        master_update_topologies();
    }
    dispatch_bubbles(MPI);
}

static const char MeshName[] = "mesh";

void Simulation:: get_meta_data( visit_handle &md ) const
{
    //! append the mesh
    {
        visit_handle mmd = _visit::mesh_meta_data(mesh, MeshName, par_size);
        VisIt_SimulationMetaData_addMesh(md, mmd);
    }
    
#if 1
    //! append the bubbles ?
    {
        const Bubble *p = bubbles.first();
        while( p )
        {
            visit_handle cmd = VISIT_INVALID_HANDLE;
            if( VisIt_CurveMetaData_alloc(&cmd) == VISIT_OKAY )
            {
                const string bubble_name = "bubble";
                VisIt_CurveMetaData_setName(cmd, bubble_name.c_str());
                
                VisIt_SimulationMetaData_addCurve(md, cmd);
            }
            p=p->next;
            break;
        }
    }
#endif
    
    //! append P on mesh
    {
        visit_handle vmd = variable_meta_data<Real>("P", MeshName);                
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }

    //! append U on mesh
    {
        visit_handle vmd = variable_meta_data<V2D>("U", MeshName);                
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
    
}

visit_handle Simulation:: get_variable( int domain, const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    
    if( name == "P" )
    {
        const int nComponents= 1;
        const int nTuples    = P.items;
        //MPI.Printf0( stderr, "Sending P: %dx%d\n", nComponents, nTuples);
        assert(P.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, P.entry);
        }
    }
    
    if( name == "U" )
    {
        const int nComponents= 2;
        const int nTuples    = U.items;
        //MPI.Printf0( stderr, "Sending U: %dx%d\n", nComponents, nTuples);
        assert(U.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, (double*)(U.entry));
        }
    }


    return h;
}


visit_handle Simulation:: get_mesh( int domain, const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if( name == MeshName )
    {
        if( VisIt_RectilinearMesh_alloc(&h) == VISIT_OKAY )
        {
            assert( h != VISIT_INVALID_HANDLE );
            const Array1D &X  = mesh.X();
            const Array1D &Y  = mesh.Y();
            
            visit_handle   hx = VISIT_INVALID_HANDLE;
            visit_handle   hy = VISIT_INVALID_HANDLE;
            VisIt_VariableData_alloc(&hx);
            VisIt_VariableData_setDataD(hx, VISIT_OWNER_SIM, 1, X.items, X.entry );
            VisIt_VariableData_alloc(&hy);
            VisIt_VariableData_setDataD(hy, VISIT_OWNER_SIM, 1, Y.items, Y.entry );
            
            VisIt_RectilinearMesh_setCoordsXY(h, hx, hy);
            
            int minRealIndex[3] = { 1,         2 ,        0 };
            int maxRealIndex[3] = { X.width-2, Y.width-3, 0 };
            
            if(  par_rank < par_size - 1)
                maxRealIndex[1]++;
            
            VisIt_RectilinearMesh_setRealIndices(h, minRealIndex, maxRealIndex);
            
        }
    }
    return h;
}

#include "yocto/ios/ocstream.hpp"

visit_handle Simulation:: get_curve( const string &name ) const
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if( name == "bubble" )
    {
        const Bubble *b = bubbles.first();
        assert(b);
        if( VisIt_CurveData_alloc(&h) == VISIT_OKAY )
        {
            // copy bubble spots coordinates
            const size_t bn = b->spots.size;
            std::cerr << "#spots=" << bn << std::endl;
            vector<Real> bx(bn,0);
            vector<Real> by(bn,0);
            const Spot *p = b->spots.head;
            
            //ios::ocstream fp("bubble.dat",false);

            for( size_t i=1; i <= bn; ++i,p=p->next)
            {
                bx[i] = p->point->vertex.x;
                by[i] = p->point->vertex.y;
                //fp("%g %g\n",bx[i],by[i]);
            }
            
            
            // make a curve
            visit_handle hcx,hcy;
            VisIt_VariableData_alloc( &hcx );
            VisIt_VariableData_alloc( &hcy );
            VisIt_VariableData_setDataD(hcx, VISIT_OWNER_COPY, 1, bn, bx());
            VisIt_VariableData_setDataD(hcy, VISIT_OWNER_COPY, 1, bn, by());
            VisIt_CurveData_setCoordsXY(h, hcx, hcy);
        
        }
    }
    return h;
}
