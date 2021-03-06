void
ui_run_clicked(void *cbdata)
{
    int prout;
    simulation_data *sim = (simulation_data *)cbdata;
    printf("ui_run_clicked from rank %d\n",rank);
   sim->runMode = SIM_RUNNING;
    VisItTimeStepChanged();
}

void
ui_dt_changed(int value, void *cbdata)
{
    simulation_data *sim = (simulation_data *)cbdata;
    printf("ui_dt_changed: %d\n", value);
}

void
ui_InitialCondition(int value, void *cbdata)
{
    simulation_data *sim = (simulation_data *)cbdata;
    printf("InitialCondition: %d\n", value);
    init_fields(value);
}

void set_interface(void *cbdata)
{
   // if(rank==0)
   // VisItUI_clicked("RUN", ui_run_clicked,cbdata);
    VisItUI_valueChanged("InitialCondition", ui_InitialCondition, cbdata);

}