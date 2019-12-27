#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "mt19937ar.h"
#include "transition_structure.h"
#include "neuronal.h"
#include "store_data.h"

double get_total_rate(int *state, state_transition *state_transitions, double *rates, int num_transitions)
{
	double total_rate = 0; 
	for(int i=0; i<num_transitions; i++){
		int state_idx = state_transitions[i].curr_state_idx;
		total_rate += state[state_idx]*rates[i];
	}
	return total_rate;
}

int get_next_state_transition(int *state, state_transition *state_transitions, double *rates, int num_transitions, double total_rate, double r)
{
	double sum = 0;
	double limit = r*total_rate;
	printf("%lf\n", total_rate);
	printf("%lf\n", r);
	for(int i=0; i<num_transitions; i++){
		int state_idx = state_transitions[i].curr_state_idx;
		sum += state[state_idx]*rates[i];
		if (sum > limit)
			return i;
	}
	return num_transitions-1;
}


void ssa(int *state, state_transition *state_transitions, double *rates, int num_transitions, double resting_v)
{
	int max_iters = 100000000;
	int iters = 0;
	double sim_time = 1000;
	double t = 0;

	init_genrand((unsigned long)time(NULL));

	double curr_v = resting_v;

	double *v_arr = (double *)calloc(max_iters, sizeof(double));
	double *t_arr = (double *)calloc(max_iters, sizeof(double));

	while(iters < max_iters){
		v_arr[iters] = curr_v;
		t_arr[iters] = t;

		double r1 = genrand_res53();
		double r2 = genrand_res53();
		double total_rate = get_total_rate(state, state_transitions, rates, num_transitions);
		
		double dt = -log(r1)/total_rate;
		int transition_idx = get_next_state_transition(state, state_transitions, rates, num_transitions,total_rate,r2);

		int curr_state_idx = state_transitions[transition_idx].curr_state_idx;
		int next_state_idx = state_transitions[transition_idx].next_state_idx;

		state[curr_state_idx]--;
		state[next_state_idx]++;

		curr_v = update_v(state,curr_v,dt);
		update_rates(rates,curr_v);	

		t = t+dt;

		printf("V : %lf\tt : %lf\n",curr_v,t);

		if(t > sim_time)
			break;

		iters++;
	}

	store_1d_arr_double("voltage_400.bin",v_arr,iters);
	store_1d_arr_double("time_400.bin",t_arr,iters);
}

int main()
{
	int num_transitions = 28;
	int num_reactants = 13;

	int *state = get_initial_state();

	double resting_v = -65;

	state_transition *neuronal_transitions = (state_transition *)calloc(num_transitions, sizeof(state_transition));

	double *rates = (double *)calloc(num_transitions, sizeof(double));

	create_neuronal_transitions(neuronal_transitions, rates, resting_v);

	for(int i=0; i<num_transitions; i++){
		printf("%d : %lf\n",i,rates[i]);
	}

	ssa(state, neuronal_transitions, rates, num_transitions, resting_v);
}


