/*****************************************************************************
 *			        Ising Model 2D			             *
 *			        Pedro H Mendes 			             *
 ****************************************************************************/

/*****************************************************************************
 *                             	   INCLUDES                                  *
 ****************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"mc.h"

/*****************************************************************************
 *                               DEFINITIONS                                 *
 ****************************************************************************/
#define 			L				20
#define 			L2 	 			(L*L)
#define 			TRAN				100000 	//1e5
#define 			TMAX				1000000 //1e6
#define 			J				1.0

/*****************************************************************************
 *                           GLOBAL VARIABLES                                *
 ****************************************************************************/
int dE, M, ET;

/*****************************************************************************
 *                              FUNCTIONS                                    *
 ****************************************************************************/
void initialize(double *boltz, int *spin, int *neigh, int *hM, int *hE, double TEMP);
void mc_routine(double *boltz, int *spin, int *neigh, int *hM, int *hE);
void print_states(int *spin, double TEMP, int* hM, int*hE, int choice);

/*****************************************************************************
 *                             MAIN PROGRAM                                  *
 ****************************************************************************/
int main(int argc, char *argv[])
{
	clock_t t_i, t_f;
	t_i = clock();

	int mcs;

	double TEMP, CPU_TIME;
	TEMP = atof(argv[1]);
	
	int *spin, *neigh;
	int *hM, *hE;
	double *boltz;
	size_t size = L2*sizeof(int); 

	spin = (int*)malloc(size);
	neigh = (int*)malloc(size);
	hM = (int*)malloc(2*size);
	hE = (int*)malloc(4*size);

	boltz = (double*)malloc(sizeof(double)*9);
	
	seed = start_randomic();

	initialize(boltz, spin, neigh, hM, hE, TEMP);

	for(mcs=0; mcs<TRAN; mcs++)
	{
		mc_routine(boltz, spin, neigh, hM, hE);
	}

#ifdef DATA	
	char Arq1[100];
	FILE *arq1;

	sprintf(Arq1, "temp_T%.3lfL%dS%ld.dat", TEMP, L, seed);
	arq1 = fopen(Arq1, "w");
	fprintf(arq1, "#seed = %ld\n#MCS\tM\tET\n", seed);
#endif

	for(mcs=0; mcs<TMAX; mcs++)
	{
		mc_routine(boltz, spin, neigh, hM, hE);
#ifdef DATA
		fprintf(arq1, "%d\t%d\t%d\n", mcs, M, ET);
#endif
	}

	print_states(spin, TEMP, hM, hE, 1);

#ifdef DATA
	fclose(arq1);
#endif

	free(spin);
	free(neigh);
	free(boltz);

	t_f = clock();
	CPU_TIME = (double)(t_f - t_i)/CLOCKS_PER_SEC;

	printf("%lf\n", CPU_TIME);

	return 0;
}

/*****************************************************************************
 *                             INITIALIZATION                                *
 ****************************************************************************/
void initialize(double *boltz, int *spin, int *neigh, int *hM, int *hE, double TEMP)
{
	int i;

	boltz[4] = exp(-4.0*J/TEMP);
	boltz[8] = exp(-8.0*J/TEMP);

	for(i=0; i<L2; i++)
	{
		if(FRANDOM < 0.5)
		{
			spin[i] = 1;
		}
		else
		{
			spin[i] = -1;
		}

		neigh[i] = 0;
	}

	ET = 0.0;
	M = 0.0;

	for(i=0; i<L2; i++)
	{
		neigh[i] += spin[(i-L+L2)%L2];
		neigh[i] += spin[(i+1)%L + (i/L)*L];
		neigh[i] += spin[(i+L)%L2];
		neigh[i] += spin[(i-1+L)%L + (i/L)*L];

		ET += spin[i]*neigh[i];
		M += spin[i];
	}

	ET = (ET*(-J))/2.0;

	for(i=0; i<2*L2; i++)
	{
		hM[i] = 0;
	}
	
	for(i=0; i<4*L2; i++)
	{
		hE[i] = 0;
	}

	return;
}

/*****************************************************************************
 *                     	      MONTE CARLO ROUTINE                            *
 ****************************************************************************/
void mc_routine(double *boltz, int *spin, int *neigh, int *hM, int *hE)
{
	int i, t;
	
	for(t=0; t<L2; t++)
	{
		i = FRANDOM*L2;

		dE = 0;
		dE = 2*spin[i]*neigh[i];

		if(dE <= 0 || FRANDOM < boltz[dE])
		{
			spin[i] *= -1;
			ET = ET + dE;
			M = M + 2*spin[i];

			neigh[(i-L+L2)%L2] += 2*spin[i];
			neigh[(i+1)%L + (i/L)*L] += 2*spin[i];
			neigh[(i+L)%L2] += 2*spin[i];
			neigh[(i-1+L)%L + (i/L)*L] += 2*spin[i];
		}
	}

	hM[M + L2] += 1;
	hE[ET + 2*L2] += 1;

	return;
}

/*****************************************************************************
 *                                PRINT STATES                               *
 ****************************************************************************/
void print_states(int *spin, double TEMP, int* hM, int*hE, int choice)
{
        int i, j;
	double mm, me;

	mm = 0.0;
	me = 0.0;

	for(i=0; i<2*L2; i++)
	{
		mm += (i-L2)*hM[i];
	}
	
	for(i=0; i<4*L2; i++)
	{
		me += (i-2*L2)*hE[i];
	}

	mm = (1.0)*mm;
	me = (1.0)*me;

	mm = mm/TMAX;
	me = me/TMAX;

        if(choice == 0)
        {
                char fp[100];
                FILE *fp1;

                sprintf(fp, "state_T%.3lfL%dS%ld.dat", TEMP, L, seed);
                fp1 = fopen(fp, "w");

                for(i=0; i<L; i++)
                {
                        for(j=0; j<L; j++)
                        {
                                fprintf(fp1, "%d ", spin[i + j*L]);
                        }
                        fprintf(fp1, "\n");
                }

                fclose(fp1);
        }
        else if(choice == 1)
        {
                char fp[100];
                FILE *fp1;

                sprintf(fp, "data_T%.3lfL%dM%.3lfE%.3lfS%ld.dat", TEMP, L, mm, me,seed);
                fp1 = fopen(fp, "w");

                for(i=0; i<L2; i++)
                {
                        fprintf(fp1, "%d ", spin[i]);
                }
                fclose(fp1);
        }

        return;
}


