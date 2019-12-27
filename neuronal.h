int *get_initial_state();
void create_neuronal_transitions(state_transition *neuronal_transitions, double *rates, double resting_v);
void update_rates(double *rates, double v);
double update_v(int *state, double v, double dt);

