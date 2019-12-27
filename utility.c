#include <stdio.h>

double get_sum(double *arr, int sz)
{
	double sum = 0;
	for(int i=0; i<sz; i++){
		sum += arr[i];
	}
	return sum;
}
