

/***************************************************************************
 * write mesh and values for all ranks
 ****************************************************************************/
void writeDomain2D(int count)
{
    int prout;
    
    DBfile *dbfile = NULL;
    
    int nz=zmax-zmin+2;
    int nx=Nx;
    int ny=1;
    int dims[] = {nx, nz};
    int ndims = 2;
    int nnodes=nx*ny*nz;
    int zf=zmax+NG;
    float *data[NC];
    
    float x[nx],z[nz];
    float *coords[] = {x, z};
    indx_t  i,j,k,p,l;
    char dirname[100];
    int cycle=count;
    double Time=count*dt;
    /* char *varName[] = {"nodal_comp0","nodal_comp1"};*/
    DBoptlist *optlist = DBMakeOptlist(2);
    DBAddOption(optlist, DBOPT_CYCLE, (void *)&cycle);
    DBAddOption(optlist, DBOPT_DTIME, (void *)&Time);
    
    

    /* prepare a rectilinear mesh. */
    for(i=0;i<nx;i++)
        x[i]=i*dx;

    for(i=0;i<nz;i++)
        z[i]=(zmin+i)*dz;
    

    sprintf(dirname, "results/bubble%d.%d", rank,count);
    dbfile = DBCreate(dirname, DB_CLOBBER, DB_LOCAL,"domain data", DB_HDF5);    
    
    
    /* write the local mesh*/
    DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, ndims,
                  DB_FLOAT, DB_COLLINEAR, optlist);
    
    
    /* Write a node-centered var. */
    for(i=0;i<(indx_t) NC;i++)
    {
        DBPutQuadvar1(dbfile, cpntName[i], "quadmesh", &fields[i][zmin][0][xmin], dims,ndims, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
    }
    

    DBClose(dbfile);
    DBFreeOptlist(optlist);
    
    
    
    /*   DBSetDir(dbfile, "..");*/
}


void writeDomain(int count)
{
    int prout;
    
    DBfile *dbfile = NULL;
    
    int nz=zmax-zmin+2;
    int nx=Nx;
    int ny=Ny;
    int dims[] = {nx, ny, nz};
    int ndims = 3;
    int nnodes=nx*ny*nz;
    int zf=zmax+NG;
    float *data[NC];
    
    float x[nx],y[ny],z[nz];
    float *coords[] = {x, y, z};
    indx_t  i,j,k,p,l;
    char dirname[100];
    int cycle=count;
    double Time=count*dt;
    /* char *varName[] = {"nodal_comp0","nodal_comp1"};*/
    DBoptlist *optlist = DBMakeOptlist(2);
    DBAddOption(optlist, DBOPT_CYCLE, (void *)&cycle);
    DBAddOption(optlist, DBOPT_DTIME, (void *)&Time);
    
    
    if(rank==size-1)
    {
        zf--;
        //   nz--;
        
        //   nnodes=nx*ny*nz;
    }
    /* prepare a rectilinear mesh. */
    for(i=0;i<nx;i++)
        x[i]=i*dx;
    for(i=0;i<ny;i++)
        y[i]=i*dy;
    for(i=0;i<nz;i++)
        z[i]=(zmin+i)*dz;
    
    /* prepare the values of the fields to send*/
    p=0;
    /*    for(k=zmax+1;k>=zmin;--k)*/  
    
    for(i=0;i<(indx_t) NC;i++)
    {
        data[i]=(float *)malloc(sizeof(float)*nnodes);
    }
    for(l=0;l<NC;l++)
    {
        for(k=zmin;k<=zf;k++)
        {
            
            for(j=ymax;j>=ymin;--j)
            {
                
                for(i=xmax;i>=xmin;--i)
                {
                    for(l=0;l<(indx_t) NC;l++)
                    {
                        data[l][p] = (float) fields[l][k][j][i];
                    }
                    p++;
                }
                
            }
            
        }
    }
    
    
    sprintf(dirname, "results/bubble%d.%d", rank,count);
    dbfile = DBCreate(dirname, DB_CLOBBER, DB_LOCAL,"domain data", DB_HDF5);    
    
    
    /* write the local mesh*/
    DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, ndims,
                  DB_FLOAT, DB_COLLINEAR, optlist);
    
    
    /* Write a node-centered var. */
    for(i=0;i<(indx_t) NC;i++)
    {
        DBPutQuadvar1(dbfile, cpntName[i], "quadmesh", data[i], dims,
                      ndims, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
    }
    
    /* DBPutQuadvar(dbfile,"nodal","quadmesh",NC,varName,data,dims,ndims,NULL, 0,DB_FLOAT,DB_NODECENT,NULL);
     */
    
    for(i=0;i<(indx_t) NC;i++)
        free(data[i]);
    
    DBClose(dbfile);
    DBFreeOptlist(optlist);
    
    
    
    /*   DBSetDir(dbfile, "..");*/
}

void write_multimesh(DBfile *dbfile, int cycle)
{
    char **meshnames = NULL;
    int dom, nmesh = size, *meshtypes = NULL;
    
    /* Create the list of mesh names. */
    meshnames = (char **)malloc(nmesh * sizeof(char *));
    for(dom = 0; dom < nmesh; ++dom)
    {
        char tmp[100];
        sprintf(tmp, "bubble%d.%d:quadmesh", dom,cycle);
        meshnames[dom] = strdup(tmp);
    }
    /* Create the list of mesh types. */
    meshtypes = (int *)malloc(nmesh * sizeof(int));
    for(dom = 0; dom < nmesh; ++dom)
        meshtypes[dom] = DB_QUAD_RECT;
    
    /* Write the multimesh. */
    DBPutMultimesh(dbfile, "quadmesh", nmesh, meshnames, meshtypes, NULL);
    
    /* Free the memory*/
    for(dom = 0; dom < nmesh; ++dom)
        free(meshnames[dom]);
    free(meshnames);
    free(meshtypes);
}



void
write_multivar(DBfile *dbfile,int cycle)
{
    int i;
    
    char **varnames = NULL;
    int dom, *vartypes = NULL;
    
    varnames = (char **)malloc(size * sizeof(char *));
    vartypes = (int *)malloc(size * sizeof(int));
    
    // Create the list of var types. 
    
    for(dom = 0; dom < size; ++dom)
        vartypes[dom] = DB_QUADVAR;
    
    
    // we write the component i of multiple domains
    for(i=0;i<(indx_t) NC;i++)
    {        
        //for each domain 
        for(dom = 0; dom < size; ++dom)
        {
            char tmp[256]={0};
            
            sprintf(tmp, "bubble%d.%d:%s", dom,cycle,cpntName[i]);
            varnames[dom] = strdup(tmp);
        }
        // Write the multivar. 
        
        DBPutMultivar(dbfile, cpntName[i], size, varnames, vartypes, NULL);
        
    }
    //  Free the memory
    
    for(dom = 0; dom < size; ++dom)
        free(varnames[dom]);
    free(varnames);
    free(vartypes);
    
}

void
write_master(int cycle)
{
    if( rank == 0 )
    {
        DBfile *dbfile = NULL;
        char fileName[100];
        
        sprintf(fileName,"./results/bubble%d.silo",cycle);
        
        /* Open the Silo file */
        dbfile = DBCreate(fileName, DB_CLOBBER, DB_LOCAL,
                          "Master file", DB_HDF5);    
        
        write_multimesh(dbfile,cycle);
        write_multivar(dbfile,cycle);
        
        /* Close the Silo file. */
        DBClose(dbfile);
    }
}
