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
} simulation_data;

#define SIM_STOPPED       0
#define SIM_RUNNING       1
void simulate_one_timestep(simulation_data *sim);

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
                    /*
                    VisItSetGetMetaData(SimGetMetaData, (void*)sim);
                    VisItSetGetMesh(SimGetMesh, (void*)sim);
                    VisItSetGetCurve(SimGetCurve, (void*)sim);
                    VisItSetGetVariable(SimGetVariable, (void*)sim);
                    VisItSetGetDomainList(SimGetDomainList, (void*)sim);
                     */
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
