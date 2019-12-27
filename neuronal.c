#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "transition_structure.h"

//Constants
#define I  0.0
#define C 1.0
#define G_L 0.3
#define G_K_MAX 36.0
#define G_NA_MAX 120.0
#define V_L -54.4
#define V_K -77.0
#define V_NA 50.0

#define K_OPEN_IDX 4
#define NA_OPEN_IDX 8

#define N_K 7200.0
#define N_NA 24000.0

int *get_initial_state()
{
	int num_reactants = 13;
	int *state = (int *)calloc(num_reactants, sizeof(int));
	
	//0 -> n0. Corresponds to the state where none of the potassium channel gates are open.
	state[0] = (int)N_K;
	
	//1 -> n1. Corresponds to the state where one of the potassium channel gates are open.
	state[1] = 0;

	//2 -> n2. Corresponds to the state where two of the potassium channel gates are open.
	state[2] = 0;

	//3 -> n3. Corresponds to the state where three of the potassium channel gates are open.
	state[3] = 0;

	//4 -> n4. Corresponds to the state where all the potassium channel gates are open.
	state[4] = 0;

	//5 -> m0h1. Corresponds to the state only the inactivating sodium channel gate is open.
	state[5] = 0;

	//6 -> m1h1. Corresponds to the state where one inactivating and one activating sodium channel gates are open.
	state[6] = 0;

	//7 -> m2h1. Corresponds to the state where one inactivating and two activating sodium channel gates are open.
	state[7] = 0;

	//8 -> m3h1. Corresponds to the state where all the sodium channel gates are open.
	state[8] = 0;

	//9 -> m0h0. Corresponds to the state where none of the sodium channel gates are open.
	state[9] = (int)N_NA;

	//10 -> m1h0. Corresponds to the state where one activating sodium channel gate is open.
	state[10] = 0;

	//11 -> m2h0. Corresponds to the state where two activating sodium channel gates are open.
	state[11] = 0;

 	//12 -> m3h0. Corresponds to the state where three activating sodium channel gates are open.
	state[12] = 0;
	return state;
}

double alpha_n(double v)
{
	double c = (v+55)/10.0;
	return (0.01*(v+55))/(1-exp(-c));
}

double alpha_m(double v)
{
	double c = (v+40)/10.0;
	return (0.1*(v+40))/(1-exp(-c));
}

double alpha_h(double v)
{
	double c = (v+65)/20.0;
	return (0.07*exp(-c));
}

double beta_n(double v)
{
	double c = (v+65)/80.0;
	return (0.125*exp(-c));
}

double beta_m(double v)
{
	double c = (v+65)/18.0;
	return (4*exp(-c));
}

double beta_h(double v)
{
	double c = (v+35)/10.0;
	return (1/(1+exp(-c)));
}

void make_transition(state_transition *state_transitions, int transition_idx, int curr_state_idx, int next_state_idx)
{
	state_transitions[transition_idx].curr_state_idx = curr_state_idx;
	state_transitions[transition_idx].next_state_idx = next_state_idx;
}

void create_neuronal_transitions(state_transition *neuronal_transitions, double *rates, double resting_v)
{

	//Transition 0 : n0 -> n1
	make_transition(neuronal_transitions, 0, 0, 1);
	rates[0] = 4*alpha_n(resting_v);

	//Transition 1 : n1 -> n0
	make_transition(neuronal_transitions, 1, 1, 0);
	rates[1] = beta_n(resting_v);

	//Transition 2 : n1 -> n2
	make_transition(neuronal_transitions, 2, 1, 2);
	rates[2] = 3*alpha_n(resting_v);

	//Transition 3 : n2 -> n1
	make_transition(neuronal_transitions, 3, 2, 1);
	rates[3] = 2*beta_n(resting_v);
	
	//Transition 4 : n2 -> n3
	make_transition(neuronal_transitions, 4, 2, 3);
	rates[4] = 2*alpha_n(resting_v);

	//Transition 5 : n3 -> n2
	make_transition(neuronal_transitions, 5, 3, 2);
	rates[5] = 3*beta_n(resting_v);

	//Transition 6 : n3 -> n4
	make_transition(neuronal_transitions, 6, 3, 4);
	rates[6] = alpha_n(resting_v);

	//Transition 7 : n4 -> n3
	make_transition(neuronal_transitions, 7, 4, 3);
	rates[7] = 4*beta_n(resting_v);

	//Transition 8 : m0h1 -> m1h1
	make_transition(neuronal_transitions, 8, 5, 6);
	rates[8] = 3*alpha_m(resting_v);

	//Transition 9 : m1h1 -> m0h1
	make_transition(neuronal_transitions, 9, 6, 5);
	rates[9] = beta_m(resting_v);

	//Transition 10 : m1h1 -> m2h1
	make_transition(neuronal_transitions, 10, 6, 7);
	rates[10] = 2*alpha_m(resting_v);

	//Transition 11 : m2h1 -> m1h1
	make_transition(neuronal_transitions, 11, 7, 6);
	rates[11] = 2*beta_m(resting_v);

	//Transition 12 : m2h1 -> m3h1
	make_transition(neuronal_transitions, 12, 7, 8);
	rates[12] = alpha_m(resting_v);

	//Transition 13 : m3h1 -> m2h1
	make_transition(neuronal_transitions, 13, 8, 7);
	rates[13] = 3*beta_m(resting_v);

	//Transition 14 : m0h0 -> m1h0
	make_transition(neuronal_transitions, 14, 9, 10);
	rates[14] = 3*alpha_m(resting_v);

	//Transition 15 : m1h0 -> m0h0
	make_transition(neuronal_transitions, 15, 10, 9);
	rates[15] = beta_m(resting_v);

	//Transition 16 : m1h0 -> m2h0
	make_transition(neuronal_transitions, 16, 10, 11);
	rates[16] = 2*alpha_m(resting_v);

	//Transition 17 : m2h0 -> m1h0
	make_transition(neuronal_transitions, 17, 11, 10);
	rates[17] = 2*beta_m(resting_v);

	//Transition 18 : m2h0 -> m3h0
	make_transition(neuronal_transitions, 18, 11, 12);
	rates[18] = 3*alpha_m(resting_v);

	//Transition 19 : m3h0 -> m2h0
	make_transition(neuronal_transitions, 19, 12, 11);
	rates[19] = beta_m(resting_v);

	//Transition 20 : m0h0 -> m0h1
	make_transition(neuronal_transitions, 20, 9, 5);
	rates[20] = alpha_h(resting_v);

	//Transition 21 : m0h1 -> m0h0
	make_transition(neuronal_transitions, 21, 5, 9);
	rates[21] = beta_h(resting_v);

	//Transition 22 : m1h0 -> m1h1
	make_transition(neuronal_transitions, 22, 10, 6);
	rates[22] = alpha_h(resting_v);

	//Transition 23 : m1h1 -> m1h0
	make_transition(neuronal_transitions, 23, 6, 10);
	rates[23] = beta_h(resting_v);

	//Transition 24 : m2h0 -> m2h1
	make_transition(neuronal_transitions, 24, 11, 7);
	rates[24] = alpha_h(resting_v);

	//Transition 25 : m2h1 -> m2h0
	make_transition(neuronal_transitions, 25, 7, 11);
	rates[25] = beta_h(resting_v);

	//Transition 26 : m3h0 -> m3h1
	make_transition(neuronal_transitions, 26, 12, 8);
	rates[26] = alpha_h(resting_v);

	//Transition 27 : m3h1 -> m3h0
	make_transition(neuronal_transitions, 27, 8, 12);
	rates[27] = beta_h(resting_v);

}


void update_rates(double *rates, double v)
{
	//Transition 0 : n0 -> n1
	rates[0] = 4*alpha_n(v);

	//Transition 1 : n1 -> n0
	rates[1] = beta_n(v);

	//Transition 2 : n1 -> n2
	rates[2] = 3*alpha_n(v);

	//Transition 3 : n2 -> n1
	rates[3] = 2*beta_n(v);
	
	//Transition 4 : n2 -> n3
	rates[4] = 2*alpha_n(v);

	//Transition 5 : n3 -> n2
	rates[5] = 3*beta_n(v);

	//Transition 6 : n3 -> n4
	rates[6] = alpha_n(v);

	//Transition 7 : n4 -> n3
	rates[7] = 4*beta_n(v);

	//Transition 8 : m0h1 -> m1h1
	rates[8] = 3*alpha_m(v);

	//Transition 9 : m1h1 -> m0h1
	rates[9] = beta_m(v);

	//Transition 10 : m1h1 -> m2h1
	rates[10] = 2*alpha_m(v);

	//Transition 11 : m2h1 -> m1h1
	rates[11] = 2*beta_m(v);

	//Transition 12 : m2h1 -> m3h1
	rates[12] = alpha_m(v);

	//Transition 13 : m3h1 -> m2h1
	rates[13] = 3*beta_m(v);

	//Transition 14 : m0h0 -> m1h0
	rates[14] = 3*alpha_m(v);

	//Transition 15 : m1h0 -> m0h0
	rates[15] = beta_m(v);

	//Transition 16 : m1h0 -> m2h0
	rates[16] = 2*alpha_m(v);

	//Transition 17 : m2h0 -> m1h0
	rates[17] = 2*beta_m(v);

	//Transition 18 : m2h0 -> m3h0
	rates[18] = 3*alpha_m(v);

	//Transition 19 : m3h0 -> m2h0
	rates[19] = beta_m(v);

	//Transition 20 : m0h0 -> m0h1
	rates[20] = alpha_h(v);

	//Transition 21 : m0h1 -> m0h0
	rates[21] = beta_h(v);

	//Transition 22 : m1h0 -> m1h1
	rates[22] = alpha_h(v);

	//Transition 23 : m1h1 -> m1h0
	rates[23] = beta_h(v);

	//Transition 24 : m2h0 -> m2h1
	rates[24] = alpha_h(v);

	//Transition 25 : m2h1 -> m2h0
	rates[25] = beta_h(v);

	//Transition 26 : m3h0 -> m3h1
	rates[26] = alpha_h(v);

	//Transition 27 : m3h1 -> m3h0
	rates[27] = beta_h(v);
}

double update_v(int *state, double v, double dt)
{
	double g_k = (state[K_OPEN_IDX]/N_K)*G_K_MAX;
	double g_na = (state[NA_OPEN_IDX]/N_NA)*G_NA_MAX;
	double dvdt = (-1/C)*(G_L*(v-V_L) + g_k*(v-V_K) + g_na*(v-V_NA) - I);

	v = v + dvdt*dt;

	return v;
}
