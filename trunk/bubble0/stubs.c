/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

/* DUMMY IMPLEMENTATIONS */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <VisItControlInterface_V2.h>
#include <VisItDataInterface_V2.h>

typedef struct
{
    int     par_rank;
    int     par_size;
    int     cycle;
    double  time;
    int     runMode;
    int     done;
    int     savingFiles;
    int     saveCounter;
    int     visitIsConnected;
} simulation_data;

#define SIM_STOPPED       0
#define SIM_RUNNING       1
void simulate_one_timestep(simulation_data *sim);

/***************************************************************************
 For interactive purposes: buttons interface for visit
 **************************************************************************/
const char *cmd_names[] = {"halt", "step", "run", "saveon","raz", "raz"}; /*For in*/
static void init_fields();
/******************************************************************************
 *
 * Purpose: This callback function returns simulation metadata.
 *
 *****************************************************************************/
visit_handle
SimGetMetaData(void *cbdata)
{
    visit_handle md = VISIT_INVALID_HANDLE;
    simulation_data *sim = (simulation_data *)cbdata;
    
    /* Create metadata. */
    if(VisIt_SimulationMetaData_alloc(&md) == VISIT_OKAY)
    {
        indx_t i;
        visit_handle mmd = VISIT_INVALID_HANDLE;
      //  visit_handle vmd = VISIT_INVALID_HANDLE;
        visit_handle vmd[NC];
        // visit_handle cmd = VISIT_INVALID_HANDLE;
        visit_handle emd = VISIT_INVALID_HANDLE;
        
        /* Set the simulation state. */
        VisIt_SimulationMetaData_setMode(md, (sim->runMode == SIM_STOPPED) ?
                                         VISIT_SIMMODE_STOPPED : VISIT_SIMMODE_RUNNING);
        VisIt_SimulationMetaData_setCycleTime(md, sim->cycle, sim->time);
        
        /* Add mesh metadata. */
        if(VisIt_MeshMetaData_alloc(&mmd) == VISIT_OKAY)
        {
            /* Set the mesh's properties.*/
            VisIt_MeshMetaData_setName(mmd, "quadmesh");
            VisIt_MeshMetaData_setMeshType(mmd, VISIT_MESHTYPE_RECTILINEAR);
           if(rmesh_dims[2]==0)
           {
            VisIt_MeshMetaData_setTopologicalDimension(mmd, 2);
            VisIt_MeshMetaData_setSpatialDimension(mmd, 2);
           }
            else
            {
                VisIt_MeshMetaData_setTopologicalDimension(mmd, 3);
                VisIt_MeshMetaData_setSpatialDimension(mmd, 3);
            }
            VisIt_MeshMetaData_setNumDomains(mmd, sim->par_size);
            VisIt_MeshMetaData_setDomainTitle(mmd, "Domains");
            VisIt_MeshMetaData_setDomainPieceName(mmd, "domain");
            VisIt_MeshMetaData_setNumGroups(mmd, 0);
            /*   VisIt_MeshMetaData_setXUnits(mmd, "cm");
             VisIt_MeshMetaData_setYUnits(mmd, "cm");
             VisIt_MeshMetaData_setZUnits(mmd, "cm");
             VisIt_MeshMetaData_setXLabel(mmd, "Width");
             VisIt_MeshMetaData_setYLabel(mmd, "Height");
             VisIt_MeshMetaData_setZLabel(mmd, "Depth");
             */
            VisIt_SimulationMetaData_addMesh(md, mmd);
            
        }
        
        for(i=0;i<NC;i++)
        {
            vmd[i]=VISIT_INVALID_HANDLE;
            // Add a variable. 
            if(VisIt_VariableMetaData_alloc(&vmd[i]) == VISIT_OKAY)
            {
                VisIt_VariableMetaData_setName(vmd[i], cpntName[i]);
                VisIt_VariableMetaData_setMeshName(vmd[i], "quadmesh");
                VisIt_VariableMetaData_setType(vmd[i], VISIT_VARTYPE_SCALAR);
                VisIt_VariableMetaData_setCentering(vmd[i], VISIT_VARCENTERING_NODE);
                VisIt_SimulationMetaData_addVariable(md, vmd[i]);
            }
            
        }
        /* Add an expression. 
         if(VisIt_ExpressionMetaData_alloc(&emd) == VISIT_OKAY)
         {
         VisIt_ExpressionMetaData_setName(emd, "zvec");
         VisIt_ExpressionMetaData_setDefinition(emd, "{zonal, zonal}");
         VisIt_ExpressionMetaData_setType(emd, VISIT_VARTYPE_VECTOR);
         
         VisIt_SimulationMetaData_addExpression(md, emd);
         }
         */
        /* Add some commands. */
        for(i = 0; i < (indx_t) (sizeof(cmd_names)/sizeof(const char *)); ++i)
        {
            visit_handle cmd = VISIT_INVALID_HANDLE;
            if(VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY)
            {
                VisIt_CommandMetaData_setName(cmd, cmd_names[i]);
                VisIt_SimulationMetaData_addGenericCommand(md, cmd);
            }
        }
    }
    
    return md;
}

visit_handle
SimGetMesh(int domain, const char *name, void *cbdata)
{
    visit_handle res = VISIT_INVALID_HANDLE;
    
    if(strcmp(name, "quadmesh") == 0)
    {
        if(VisIt_RectilinearMesh_alloc(&res) != VISIT_ERROR)
        {
            if(rmesh_dims[2]==0) // we are in 2D
            {
                int i,minRealIndex[2]={0,0}, maxRealIndex[2]={0,0};
                float *rmesh[2];
                visit_handle h[2];
                 fprintf(stderr,"proc %d\t:simgetmesh: %d\t%d\t%d in 2D\n",rank,rmesh_dims[0],rmesh_dims[1],rmesh_dims[2]);
                 fflush(stderr);
                
                for(i=0;i<2;i++)
                {
                    minRealIndex[i] = 0;
                    maxRealIndex[i] = rmesh_dims[i]-1;
                    // attention do not free ! since we set the flag VISIT_OWNER_VISIT
                    rmesh[i] = (float *)malloc(sizeof(float) * rmesh_dims[i]);
                }
                
                for(i = 0; i < rmesh_dims[0]; ++i)
                    rmesh[0][i] = i*dx;
                for(i = 0; i < rmesh_dims[1]; ++i)
                   rmesh[1][i] = (zmin+i)*dz;
                
                for(i=0;i<2;i++)
                {
                    VisIt_VariableData_alloc(&h[i]);
                    VisIt_VariableData_setDataF(h[i], VISIT_OWNER_VISIT, 1,rmesh_dims[i], rmesh[i]);
                    
                }
                VisIt_RectilinearMesh_setCoordsXY(res, h[0], h[1]);
                VisIt_RectilinearMesh_setRealIndices(res, minRealIndex, maxRealIndex);
            }
            else
            {
                int i,minRealIndex[3]={0,0,0}, maxRealIndex[3]={0,0,0};
                float *rmesh[3];
                visit_handle h[3];
                
                // fprintf(stderr,"proc %d\t:simgetmesh: %d\t%d\t%d\n",rank,rmesh_dims[0],rmesh_dims[1],rmesh_dims[2]);
                // fflush(stderr);
                
                for(i=0;i<3;i++)
                {
                    minRealIndex[i] = 0;
                    maxRealIndex[i] = rmesh_dims[i]-1;
                    // attention do not free !
                    rmesh[i] = (float *)malloc(sizeof(float) * rmesh_dims[i]);
                }
                
                for(i = 0; i < rmesh_dims[0]; ++i)
                    rmesh[0][i] = i*dx;
                for(i = 0; i < rmesh_dims[1]; ++i)
                    rmesh[1][i] = i*dy;
                for(i = 0; i < rmesh_dims[2]; ++i)
                    rmesh[2][i] =(zmin+i)*dz;
                
                
                for(i=0;i<3;i++)
                {
                    VisIt_VariableData_alloc(&h[i]);
                    VisIt_VariableData_setDataF(h[i], VISIT_OWNER_VISIT, 1,rmesh_dims[i], rmesh[i]);
                    
                }
                VisIt_RectilinearMesh_setCoordsXYZ(res, h[0], h[1],h[2]);
                VisIt_RectilinearMesh_setRealIndices(res, minRealIndex, maxRealIndex);
            }
                
             
        }
        else
        {
            fprintf(stderr,"proc: %d:Erreur allocation in SimGetMesh\n",rank);
        }
    }
    
    return res;
}
/*
visit_handle
SimGetVariableWorking(int domain, const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    //  simulation_data *sim = (simulation_data *)cbdata;
    
    // fprintf(stderr,"proc %d: SimGetVariable\n",rank);
    
    if(strcmp(name, "u") == 0)
    {
        float *zoneptr;
        float  *rmesh_zonal;
        int i, j, k, nTuples;
        
        
        // Calculate a zonal variable that moves around. 
        rmesh_zonal = (float*)malloc(sizeof(float) * (rmesh_dims[0]-1) * (rmesh_dims[1]-1)*(rmesh_dims[2]-1));
        zoneptr = rmesh_zonal;
        
        for(k=zmin;k<zmax;k++)
        {
            for(j = 0; j < rmesh_dims[1]-1; ++j)
            {
                for(i = 0; i < rmesh_dims[0]-1; ++i)
                {
                    
                    *zoneptr++ = fields[0][k][j][i];
                }
            }
        }
        nTuples = (rmesh_dims[0]-1) * (rmesh_dims[1]-1)*(rmesh_dims[2]-1);
        VisIt_VariableData_alloc(&h);
        VisIt_VariableData_setDataF(h, VISIT_OWNER_VISIT, 1,
                                    nTuples, rmesh_zonal);
    }
    
    return h;
}
*/
visit_handle
SimGetVariable(int domain, const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    int toPlot=0;
    int i;
    //  simulation_data *sim = (simulation_data *)cbdata;
    
    // fprintf(stderr,"proc %d: SimGetVariable\n",rank);
    for(i=0;i<NC;i++)
       if(strcmp(name, "u") == 0)
           toPlot=i;
        
   // if(strcmp(name, "u") == 0)
     if(rmesh_dims[2]==0) // we are in 2D
    {
        float *zoneptr;
        float  *rmesh_zonal;
        int i, j, k, nTuples;
        if(rmesh_dims[2]==0) //we are in 2D
            nTuples = (rmesh_dims[0]) * (rmesh_dims[1]);

        else
            nTuples = (rmesh_dims[0]) * (rmesh_dims[1])*(rmesh_dims[2]);

        
        // Calculate a zonal variable that moves around. 
        rmesh_zonal = (float*)malloc(sizeof(float)*nTuples);
        zoneptr = rmesh_zonal;
        
        // A RETRAVAILLER
        if((size==2)&&rank==0)      
        {
           // for(k=zmin;k<=zmax;k++)
                for(k=zmax+NG;k>=zmin;--k)
            {
                j=0;
                // for(j = 0; j < rmesh_dims[1]; ++j)
                {
                    for(i = 0; i < rmesh_dims[0]; ++i)
                    {
                        
                        *zoneptr++ = fields[toPlot][k][j][i];
                    }
                }
            }
        }
        else
        {
            for(k=zmin;k<=zmax+NG;k++)
            {
                j=0;
                // for(j = 0; j < rmesh_dims[1]; ++j)
                {
                    for(i = 0; i < rmesh_dims[0]; ++i)
                    {
                        
                        *zoneptr++ = fields[toPlot][k][j][i];
                    }
                }
            }
        }
        VisIt_VariableData_alloc(&h);
        VisIt_VariableData_setDataF(h,VISIT_OWNER_VISIT,1,nTuples, rmesh_zonal);
    }
    else // we are in 3D
    {
        float *zoneptr;
        float  *rmesh_zonal;
        int i, j, k, nTuples;
        nTuples = (rmesh_dims[0]) * (rmesh_dims[1])*(rmesh_dims[2]);
        
        
        // Calculate a zonal variable that moves around. 
        rmesh_zonal = (float*)malloc(sizeof(float)*nTuples);
        zoneptr = rmesh_zonal;
        
        // A RETRAVAILLER
        if((size==2)&&rank==0)      
        {
            // for(k=zmin;k<=zmax;k++)
            for(k=zmax+NG;k>=zmin;--k)
            {
                for(j = 0; j < rmesh_dims[1]; ++j)
                {
                    for(i = 0; i < rmesh_dims[0]; ++i)
                    {
                        
                        *zoneptr++ = fields[toPlot][k][j][i];
                    }
                }
            }
        }
        else
        {
            for(k=zmin;k<=zmax+NG;k++)
            {
                for(j = 0; j < rmesh_dims[1]; ++j)
                {
                    for(i = 0; i < rmesh_dims[0]; ++i)
                    {
                        
                        *zoneptr++ = fields[toPlot][k][j][i];
                    }
                }
            }
        }
        VisIt_VariableData_alloc(&h);
        VisIt_VariableData_setDataF(h,VISIT_OWNER_VISIT,1,nTuples, rmesh_zonal);
    }
    
    return h;
}

visit_handle
SimGetDomainList(const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    
    //fprintf(stderr,"proc %d: SimGetDomainList\n",rank);
    
    if(VisIt_DomainList_alloc(&h) != VISIT_ERROR)
    {
        visit_handle hdl;
        int *iptr = NULL;
        simulation_data *sim = (simulation_data *)cbdata;
        
        iptr = (int *)malloc(sizeof(int));
        *iptr = sim->par_rank;
        
        VisIt_VariableData_alloc(&hdl);
        VisIt_VariableData_setDataI(hdl, VISIT_OWNER_VISIT, 1, 1, iptr);
        VisIt_DomainList_setDomains(h, sim->par_size, hdl);
    }
    return h;
}
void
simulation_data_ctor(simulation_data *sim,int rank,int size)
{
    sim->par_rank = rank;
    sim->par_size = size;
    sim->cycle = 0;
    sim->time = 0.;
    sim->runMode = SIM_STOPPED;
    sim->done = 0;
    sim->savingFiles = 0;
    sim->saveCounter = 0;
    sim->visitIsConnected=0;
}




static int visit_broadcast_int_callback(int *value, int sender)
{
    return MPI_Bcast(value, 1, MPI_INT, sender, MPI_COMM_WORLD);
}

static int visit_broadcast_string_callback(char *str, int len, int sender)
{
    return MPI_Bcast(str, len, MPI_CHAR, sender, MPI_COMM_WORLD);
}

void
simulation_data_dtor(simulation_data *sim)
{
}


/* CHANGE 1 */
#define VISIT_COMMAND_PROCESS 0
#define VISIT_COMMAND_SUCCESS 1
#define VISIT_COMMAND_FAILURE 2

/* Helper function for ProcessVisItCommand */
static void BroadcastSlaveCommand(int *command)
{
    MPI_Bcast(command, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

/* Callback involved in command communication. */
void SlaveProcessCallback()
{
    int command = VISIT_COMMAND_PROCESS;
    BroadcastSlaveCommand(&command);
}

void read_input_deck(simulation_data *sim)
{
    /* Read in problem setup. */
}

void write_vis_dump(simulation_data *sim)
{
    /* Write visualization dump. */
}
/******************************************************************************
 * Purpose: Process commands from viewer on all processors. 
 *****************************************************************************/
int ProcessVisItCommand(simulation_data *sim)
{
    int command;
    if (sim->par_rank==0)
    {  
        int success = VisItProcessEngineCommand();
        
        if (success)
        {
            command = VISIT_COMMAND_SUCCESS;
            BroadcastSlaveCommand(&command);
            return 1;
        }
        else
        {
            command = VISIT_COMMAND_FAILURE;
            BroadcastSlaveCommand(&command);
            return 0;
        }
    }
    else
    {
        /* Note: only through the SlaveProcessCallback callback
         * above can the rank 0 process send a VISIT_COMMAND_PROCESS
         * instruction to the non-rank 0 processes. */
        while (1)
        {
            BroadcastSlaveCommand(&command);
            switch (command)
            {
                case VISIT_COMMAND_PROCESS:
                    VisItProcessEngineCommand();
                    break;
                case VISIT_COMMAND_SUCCESS:
                    return 1;
                case VISIT_COMMAND_FAILURE:
                    return 0;
            }
        }
    }
    return 1;
}

/**************************************************************************
 *Callback function for control commands, which are the buttons in the 
 * GUI's Simulation window. This type of command is handled automatically
 * provided that you have registered a command callback such as this.
**************************************************************************/
void ControlCommandCallback(const char *cmd, const char *args, void *cbdata)
{
    simulation_data *sim = (simulation_data *)cbdata;
    
    if(strcmp(cmd, "halt") == 0)
        sim->runMode = SIM_STOPPED;
    else if(strcmp(cmd, "step") == 0)
        simulate_one_timestep(sim);
    else if(strcmp(cmd, "run") == 0)
        sim->runMode = SIM_RUNNING;
    else if(strcmp(cmd, "raz") == 0)
        init_fields();
    else if(strcmp(cmd, "addplot") == 0)
    {
        VisItExecuteCommand("AddPlot(\"Pseudocolor\", \"zonal\")\n");
        VisItExecuteCommand("DrawPlots()\n");
    }
}

/* Called to handle case 3 from VisItDetectInput where we have console
 * input that needs to be processed in order to accomplish an action.
 */
void
ProcessConsoleCommand(simulation_data *sim)
{
    /* Read A Command */
    char cmd[1000];
    
    if (sim->par_rank == 0)
    {
        int iseof = (fgets(cmd, 1000, stdin) == NULL);
        if (iseof)
        {
            sprintf(cmd, "quit");
            printf("quit\n");
        }
        
        if (strlen(cmd)>0 && cmd[strlen(cmd)-1] == '\n')
            cmd[strlen(cmd)-1] = '\0';
    }
    
    /* Broadcast the command to all processors. */
    MPI_Bcast(cmd, 1000, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    if(strcmp(cmd, "quit") == 0)
        sim->done = 1;
    else if(strcmp(cmd, "halt") == 0)
        sim->runMode = SIM_STOPPED;
    else if(strcmp(cmd, "step") == 0)
        simulate_one_timestep(sim);
    else if(strcmp(cmd, "run") == 0)
        sim->runMode = SIM_RUNNING;
    else if(strcmp(cmd, "raz") == 0)
        init_fields();
    else if(strcmp(cmd, "update") == 0)
    {
        VisItTimeStepChanged();
        VisItUpdatePlots();
    }
    else if(strcmp(cmd, "saveon") == 0)
        sim->savingFiles = 1;
    else if(strcmp(cmd, "saveoff") == 0)
        sim->savingFiles = 0;
    else if(strcmp(cmd, "addplot") == 0)
    {
        VisItExecuteCommand("AddPlot(\"Pseudocolor\", \"zonal\")\n");
        VisItExecuteCommand("DrawPlots()\n");
    }
}

void mainloop(simulation_data *sim)
{
    int blocking, visitstate, err = 0;
    
    /* If we're not running by default then simulate once there's something
     * once VisIt connects.
     */
    if(sim->runMode == SIM_STOPPED)
        simulate_one_timestep(sim);
    
    if (sim->par_rank == 0)
    {
        fprintf(stderr, "command> ");
        fflush(stderr);
    }
    
    do
    {
        blocking = (sim->runMode == SIM_RUNNING) ? 0 : 1;
        /* Get input from VisIt or timeout so the simulation can run. */
        if(sim->par_rank == 0)
        {
            visitstate = VisItDetectInput(blocking, fileno(stdin));
        }
        /* Broadcast the return value of VisItDetectInput to all procs. */
        MPI_Bcast(&visitstate, 1, MPI_INT, 0, MPI_COMM_WORLD);
        /* Do different things depending on the output from VisItDetectInput. */

        switch(visitstate)
        {
            case 0:
                /* There was no input from VisIt, return control to sim. */
                simulate_one_timestep(sim);
                break;
            case 1:
                /* VisIt is trying to connect to sim. */
                if(VisItAttemptToCompleteConnection() == VISIT_OKAY)
                {
                    fprintf(stderr, "VisIt connected\n");
                    VisItSetCommandCallback(ControlCommandCallback, (void*)sim);
                    VisItSetSlaveProcessCallback(SlaveProcessCallback);
                    
                    VisItSetGetMetaData(SimGetMetaData, (void*)sim);
                    VisItSetGetMesh(SimGetMesh, (void*)sim);
                   //VisItSetGetCurve(SimGetCurve, (void*)sim);
                    VisItSetGetVariable(SimGetVariable, (void*)sim);
                    VisItSetGetDomainList(SimGetDomainList, (void*)sim);
                    sim->visitIsConnected=1;
                     
                }
                else 
                {
                    /* Print the error message */
                    char *err = VisItGetLastError();
                    fprintf(stderr, "VisIt did not connect: %s\n", err);
                    free(err);
                }
                break;
            case 2:
                /* VisIt wants to tell the engine something. */
                if(!ProcessVisItCommand(sim))
                {
                    /* Disconnect on an error or closed connection. */
                    VisItDisconnect();
                    sim->visitIsConnected=0;

                    /* Start running again if VisIt closes. */
                    /*sim->runMode = SIM_RUNNING;*/
                }
                break;
            case 3:
                /* VisItDetectInput detected console input - do something with it.
                 * NOTE: you can't get here unless you pass a file descriptor to
                 * VisItDetectInput instead of -1.
                 */
                ProcessConsoleCommand(sim);
                if (sim->par_rank == 0)
                {
                    fprintf(stderr, "command> ");
                    fflush(stderr);
                }
                break;
            default:
                fprintf(stderr, "Can't recover from error %d!\n", visitstate);
                err = 1;
                break;
        }
    } while(!sim->done && err == 0);
}
