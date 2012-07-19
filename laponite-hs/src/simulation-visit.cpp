#include "simulation.hpp"

static const char MeshName[] = "mesh";

void Simulation:: get_meta_data( visit_handle &md ) const
{
    
    //! append an UI command
    {
        visit_handle cmd = VISIT_INVALID_HANDLE;
        if(VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY)
        {
            VisIt_CommandMetaData_setName(cmd, "save1" );
            VisIt_SimulationMetaData_addGenericCommand(md, cmd);
        }
    }
    
    //! append an UI command
    {
        visit_handle cmd = VISIT_INVALID_HANDLE;
        if(VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY)
        {
            VisIt_CommandMetaData_setName(cmd, "ins" );
            VisIt_SimulationMetaData_addGenericCommand(md, cmd);
        }
    }
    
    //! append the mesh
    {
        visit_handle mmd = _visit::mesh_meta_data(mesh, MeshName, par_size);
        VisIt_SimulationMetaData_addMesh(md, mmd);
    }
    
    if( !parallel )    //! append the bubbles ?
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
    
    //! append P on mesh
    {
        visit_handle vmd = variable_meta_data<Real>("P", MeshName);                
        VisIt_SimulationMetaData_addVariable(md, vmd);
    }
    
    //! append B on mesh
    {
        visit_handle vmd = variable_meta_data<Real>("B", MeshName);                
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
        return h;
    }
    
    if( name == "B" )
    {
        const int nComponents= 1;
        const int nTuples    = B.items;
        //MPI.Printf0( stderr, "Sending B: %dx%d\n", nComponents, nTuples);
        assert(B.entry!=NULL);
        if(VisIt_VariableData_alloc(&h) == VISIT_OKAY)
        {
            VisIt_VariableData_setDataD(h, VISIT_OWNER_SIM, nComponents, nTuples, B.entry);
        }
        return h;
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
        return h;
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
            //const Array1D &X  = mesh.X();
            //const Array1D &Y  = mesh.Y();
            
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
            const size_t bn = b->size+1;
            vector<Real> bx(bn,0);
            vector<Real> by(bn,0);
            const Point *p = b->root;
            
            //ios::ocstream fp("bubble.dat",false);
            
            for( size_t i=1; i <= bn; ++i,p=p->next)
            {
                bx[i] = p->vertex.x;
                by[i] = p->vertex.y;
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


void Simulation:: perform(const string &cmd)
{
    if( cmd == "save1" )
    {
        if( master )
        {
            bubbles.first()->save_vtk("bubble.vtk");
            bubbles.first()->save_vtk_t("bubble_t.vtk");
            bubbles.first()->save_vtk_n("bubble_n.vtk");
        }
    }
    
    if( cmd == "ins" )
    {
        if( master )
        {
            vector<V2D> pts;
            collect_inside(pts);
            bubbles.first()->save_dat("bubble.dat");
            ios::ocstream fp( "inside.dat", false );
            for( size_t i=pts.size(); i>0; --i )
            {
                fp("%g %g\n",pts[i].x, pts[i].y);
            }
        }
    }
}

void Simulation:: step()
{
    VisIt::Simulation::step();
    // move concerned points
    advect_points(1);
    
    // send back information to master
    assemble_bubbles(MPI);
    
    //process topologies on master
    check_and_dispatch_bubbles();
    
    
    init_exchange();
    wait_exchange();
    
}


