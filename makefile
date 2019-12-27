objects = mt19937ar.o neuronal.o store_data.o gillespie.o

neuronal_stochastic : $(objects)
	gcc -o $@ $(objects) -lm

mt19937ar.o : mt19937ar.c
	gcc -c mt19937ar.c

neuronal.o : neuronal.c
	gcc -c neuronal.c -lm

store_data.o : store_data.c
	gcc -c store_data.c

gillespie.o : mt19937ar.c neuronal.c gillespie.c
	gcc -c gillespie.c -lm
