#include <stdio.h>
#include <stdlib.h>

void store_1d_arr(char *fname, int *arr, int sz)
{
	FILE *fp = fopen(fname,"wb");
	if(fp == NULL){
		printf("Unable to create File.\n");
		exit(0);
	}
	fwrite(arr,sizeof(int),sz,fp);
	fclose(fp);
}

void store_1d_arr_double(char *fname, double *arr, int sz)
{
	FILE *fp = fopen(fname,"wb");
	if(fp == NULL){
		printf("Unable to create File.\n");
		exit(0);
	}
	fwrite(arr,sizeof(double),sz,fp);
	fclose(fp);
}

void store_2d_arr(char *fname, int **arr, int sz)
{
	FILE *fp = fopen(fname,"wb");
	if(fp == NULL){
		printf("Unable to create File.\n");
		exit(0);
	}
	for(int i = 0; i<sz; i++)
		fwrite(arr[i],sizeof(int),sz,fp);
	fclose(fp);
}

void store_2d_arr_pack(char *fname, int ***arr_pack, int sz, int pack_sz)
{
	FILE *fp = fopen(fname,"wb");
	if(fp == NULL){
		printf("Unable to create File.\n");
		exit(0);
	}
	for(int i = 0; i<pack_sz; i++){
		for(int j=0; j<sz; j++)
			fwrite(arr_pack[i][j],sizeof(int),sz,fp);
	}
	fclose(fp);
}

void store_1d_arr_dat(char *fname, int *arr, int sz)
{
	FILE *fp = fopen(fname,"w");
	if(fp == NULL){
		printf("Unable to create File.\n");
		exit(0);
	}
	for(int i=0; i<sz; i++){
		fprintf(fp,"%d\n",arr[i]);
	}
	fclose(fp);
}

