/**************************************************************************
 *	  		     Potts Model 2D				  *
 *************************************************************************/

/**************************************************************************
 *	 		       INCLUDES 				  *
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mc.h"

/**************************************************************************
 *			      DEFINITIONS				  *
 *************************************************************************/
#define 			L			64	
#define 			L2			(L*L)
#define 			J		 	1.
#define 			KB			1.
#define 			TRAN			100000
#define 			TMAX			1000000

/**************************************************************************
 *	 		   GLOBAL VARIABLES				  *
 *************************************************************************/
int M, ET;
double exp_factor;

/**************************************************************************
 *  			     FUNCTIONS 					  *
 *************************************************************************/
void initialize(int *spin, int **neigh, double *prob_met, double *prob_hb, int *hE, int **hM, int _Q, double TEMP);
void sweep_heatbath(int *spin, int **neigh, double *prob_hb, int _Q);
void sweep_metropolis(int *spin, int **neigh, double *prob_met, int _Q);
void states(int *spin, int **neigh, int *hE, int **hM, int choice, int step, int _Q);
void print_states(int *spin, double TEMP, int *hE, int **hM, int _Q, int choice);
void swendsen_wang_step(int *spin, int **neigh, int _Q);
void wolff_step(int *spin, int **neigh, int _Q);
void gnuplot_view(int tempo, int *s);

/**************************************************************************
 *		  	    MAIN PROGRAM				  *
 *************************************************************************/
int main(int argc, char *argv[])
{
	int *spin, **neigh;
	int *hE, **hM;
	double *prob_met,*prob_hb;
	int i, mcs;
	int Q;	
	double TEMP;

	 // leio os parametros passados pela linha de comando
	fprintf(stderr,"# argc: %d\n",argc);
	if(argc==5)
	{
		for(i=1;i<argc;i++)
		{
			if(!strcmp(argv[i],"-Q"))
			{
				Q=atoi(argv[++i]);
			}
			else if(!strcmp(argv[i],"-T"))
			{
				TEMP=atof(argv[++i]);
			}
			else
			{
				fprintf(stderr,"Error.	Argument '%s' is not recognized.\n", argv[i]);
				exit(-1);
			}
		}
	}
	else
	{
		fprintf(stderr,"Error.	Number of Arguments is wrong!!.\n");
		exit(-1);
	}
	fprintf(stderr,"\n Simulacao do modelo de Potts com:\n");
	fprintf(stderr," Lsize	 = %d\n",L);
	fprintf(stderr," Q			 = %d\n",Q);
	fprintf(stderr," Temp		= %f\n",TEMP);
	fprintf(stderr," Tc		= %f\n\n",1./log(1.0+sqrt(Q)));
		
	// fator pro SW
	exp_factor = 1-exp(-1./TEMP);

	// alocacao de matrizes 1D
	spin = (int*)malloc(L2*sizeof(int));
	prob_met = (double*)malloc(5*sizeof(double));
	prob_hb = (double*)malloc(5*sizeof(double));
	hE = (int*)malloc(TMAX*sizeof(int));
	// alocacao de matrizes 2D
	//
	neigh = (int**)malloc(L2*sizeof(int*));
	for(i=0; i<L2; i++)
	{
		neigh[i] = (int*)malloc(4*sizeof(int));
	}
	//
	hM = (int**)malloc(TMAX*sizeof(int*));
	for(i=0; i<TMAX; i++)
	{
		hM[i] = (int*)malloc(Q*sizeof(int));
	}

	// inicializo o RNG
	seed = start_randomic();

	initialize(spin, neigh, prob_met, prob_hb, hE, hM, Q, TEMP);
	states(spin, neigh, hE, hM, 0, 0,Q);

	for(mcs=0; mcs<TRAN; mcs++)
	{
		//sweep_heatbath(spin, neigh, prob_hb, Q);i
		//sweep_metropolis(spin, neigh, prob_met, Q);	
		swendsen_wang_step(spin, neigh, Q);
		//wolff_step(spin, neigh, Q);
		states(spin, neigh, hE, hM, 0, mcs,Q);
#ifdef GNUPLOT
		gnuplot_view(mcs,spin);
#endif
	}

#ifdef DATA
	char Arq[100];
	FILE *arq;
	sprintf(Arq, "potts_Q%d_T%lf_L%d_S%ld.dsf", Q, TEMP, L, seed);
	arq = fopen(Arq, "w");
	fprintf(arq, "#SEED: %ld\n#MCS\tM\tET\n", seed);
#endif
			
	for(mcs=0; mcs<TMAX; mcs++)
	{
		//sweep_heatbath(spin, neigh, prob_hb, Q);
		//sweep_metropolis(spin, neigh, prob_met, Q);
		swendsen_wang_step(spin, neigh, Q);
		//wolff_step(spin, neigh, Q);
		states(spin, neigh, hE, hM, 1, mcs,Q);
#ifdef GNUPLOT
		gnuplot_view(mcs,spin);
#endif
#ifdef DATA
		fprintf(arq, "%d\t%d\t%d\n", mcs, M, ET);
#endif
	}

	print_states(spin, TEMP, hE, hM, Q, 0);

#ifdef DATA
	fclose(arq);
#endif

	// libero as matrizes alocadas 1D
	free(spin);
	free(prob_met);
	free(prob_hb);
	free(hE);
	// libero as matrizes alocadas 2D
	
	for(i=0; i<L2; i++)
	{
		free(neigh[i]);
	}
	free(neigh);
	
	for(i=0; i<Q; i++)
	{
		free(hM[i]);
	}
	free(hM);

	return 0;
}

/**************************************************************************
 *			   INITIALIZATION				  *
 *************************************************************************/
void initialize(int *spin, int **neigh, double *prob_met, double *prob_hb, int *hE, int **hM, int _Q, double TEMP)
{
	int i,j;

	// spins sempre aleatorios
	for(i=0; i<L2; i++)
	{
		spin[i]=FRANDOM*_Q;
	}
 	// spins sempre ordenados
	for(i=0; i<L2; i++)
	{
		spin[i]=1;
	}

	// matriz de vizinhos
	for(i=0; i<L2; i++)
	{
		neigh[i][0] = (i-L+L2)%L2;		//up
		neigh[i][1] = (i+1)%L + (i/L)*L;	//right
		neigh[i][2] = (i+L)%L2;		//down
		neigh[i][3] = (i-1+L)%L + (i/L)*L;	//left
	}

	// matrizes series temporais
	for(i=0; i<TMAX; i++)
	{
		hE[i] = 0;
		for(j=0; j<_Q; j++)
		{
			hM[i][j] = 0;
		}
	}

	// matriz de prob
	for(i=0; i<5; i++)
	{
		prob_hb[i] = exp(i/TEMP);	
		prob_met[i] = exp(-i/TEMP);	
	}

	return;
}

/**************************************************************************
 *		 	MONTE CARLO ROUTINE				  *
 *	Heat bath algorithm						  *
 *************************************************************************/
void sweep_heatbath(int *spin, int **neigh, double *prob_hb, int _Q)
{
	int site, i, j;
	double ran, norm;
	int howmany[_Q];

	// MCS
	for(i=0; i<L2; i++)
	{
		// escolho um sitio aleatorio
		site = (int)(FRANDOM*L2);
		// conto quantos dos seus vizinhos sao de cada Q
		// inicializo
		for(j=0; j<_Q; j++)
		{
			howmany[j] = 0;
		}
		// conto
		++howmany[spin[neigh[site][0]]];
		++howmany[spin[neigh[site][1]]];
		++howmany[spin[neigh[site][2]]];
		++howmany[spin[neigh[site][3]]];
		// tower sampling
		norm = 0.0;
		// acumulo as probs
		for(j=0; j<_Q; j++)
		{
			// chance maior de flipar para o estado mais presente nos vizinhos
			norm += prob_hb[howmany[j]];
		}
		// posicao aleatoria na torre
		ran = FRANDOM*norm;
		// procuro a caixa correspondente a qual Q
		for(j=0; j<_Q; j++)
		{
			if(ran < prob_hb[howmany[j]])
			{
				spin[site] = j;
				break;
			}
			else
			{
				ran -= prob_hb[howmany[j]];
			}
		}
	}
	
	return;
}

/**************************************************************************
 *			  MONTE CARLO ROUTINE				  *
 *	Metropolis algorithm						  *
 *************************************************************************/
void sweep_metropolis(int *spin, int **neigh, double *prob_met, int _Q)
{
	int site, i, j;
//	double ran, norm;
	int dQ,newspin,ssite,e0,ef,deltae;

	// MCS
	for(i=0; i<L2; i++)
	{
		// escolho um sitio aleatorio
		site = (int)(FRANDOM*L2);
		ssite = spin[site];
		// escolho um novo Q
		// acho que aleatorio funcio mas nao faz sentido
		//newspin = FRANDOM*_Q;
		dQ = (_Q-1)*FRANDOM; // supondo que nao gero o 1 !!!
		//dQ = (_Q-2)*FRANDOM;
		newspin = (spin[site] + 1 + dQ)%_Q;
		//newspin = (ssite+1)%2;
		if((newspin==_Q)||(newspin==ssite))
		{
			fprintf(stderr,"erro no newspin: %d %d (METROPOLIS)\n",newspin,dQ);
			fprintf(stderr,"erro no newspin: %d %d (METROPOLIS)\n",newspin,ssite);
			exit(-1);
		}
		// energia inicial e final
		
		e0 = 0;
		ef = 0;
		
		for(j=0;j<4;j++)
		{
			e0 -= (ssite==spin[neigh[site][j]]);
			ef -= (newspin==spin[neigh[site][j]]);
		}
		deltae = ef - e0;
		
		if((deltae<=0)||(FRANDOM<prob_met[deltae]))
		{
			spin[site] = newspin;
			//ET += deltae;
		}
	}
	
	return;
}

/**************************************************************************
 *				 STATES		 			  *
 *************************************************************************/
void states(int *spin, int **neigh, int *hE, int **hM, int choice, int step, int _Q)
{
	int i, j;

	ET=0;
	M=0;

	if(choice == 1)
	{
		for(j=0; j<_Q; j++)
		{
			hM[step][j] = 0;
		}
	}

	for(i=0; i<L2; i++)
	{
		// so o de baixo e o da direita
		for(j=1; j<3; j++)
		{
			if(spin[i] == spin[neigh[i][j]])
			{
				ET--;
			}
		}
		if(choice == 1)
		{
			hM[step][spin[i]]++;
		}
		// esse M nao esta fazendo sentido
		M += spin[i];
	}

	// ET deveria ser dividido por 2
	if(choice == 1)
	{
		hE[step] = ET;
	}
	
	return;
}

/**************************************************************************
 *				PRINT STATES				  *
 *************************************************************************/
void print_states(int *spin, double TEMP, int *hE, int **hM, int _Q, int choice)
{
	int i, j;
	double mm, me;

	mm = 0.0;
	me = 0.0;

	for(i=0; i<TMAX; i++)
	{
		mm += hM[i][0];//hM[i];
		me += hE[i];
	}

	mm = (1.0)*mm;
	me = (1.0)*me;

	mm = mm/TMAX;
	me = me/TMAX;

	if(choice == 0)
	{
		char fp[100];
		FILE *fp1;

		sprintf(fp, "state_T%lfQ%dL%dS%ld.dsf", TEMP, _Q, L, seed);
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

		sprintf(fp, "data_T%lfQ%dL%dM%.3lfE%.3lfS%ld.dsf",TEMP, _Q, L, mm, me,seed);
		fp1 = fopen(fp, "w");

		for(i=0; i<L2; i++)
		{
			fprintf(fp1, "%d ", spin[i]);
		}
		fclose(fp1);
	}

	return;
}

/**************************************************************************
 *		 	 	SWENDSEN WANG				  *
 *************************************************************************/
void swendsen_wang_step(int *spin, int **neigh, int _Q)
{
	int i, j, k;
	unsigned long *label;
	signed int *flip;
	double temp;

	label = malloc(L2*sizeof(unsigned long));
	flip = malloc(L2*sizeof(signed int));

	for(i=0; i<L2; i++)
	{
		label[i] = i;
		flip[i] = -1;
	}

	for(i=0; i<L2; i++)
	{
		if( (spin[i]==spin[neigh[i][1]]) && (FRANDOM < exp_factor) )
		{
			unionfind(i,neigh[i][1],label);
		}
		if( (spin[i]==spin[neigh[i][2]]) && (FRANDOM < exp_factor) )
		{
			unionfind(i,neigh[i][2],label);
		}
	}
	
	for(i=0; i<L2; i++)
	{
		j = i;
	
		while(label[j] != j)
		{
			j = label[j];
		}
	
		if(flip[j] == -1)
		{
			temp = FRANDOM;
			for(k=0; k<_Q; k++)
			{
				if(temp<1./_Q)
				{
					flip[j] = k;
					break;
				}
				else
				{
					temp -= 1./_Q;
				}
			}
		}
		spin[i] = flip[j];
	}
	
	free(label);
	free(flip);

	return;
}



/**************************************************************************
 *			 	    WOLFF				  *
 *	baseado no isng do barkema					  *
 *	conferir se esta certo pro potts..				  *
 *************************************************************************/
void wolff_step(int *spin, int **neigh, int _Q) {

	int i,j,sp,oldspin,newspin,current,nn,dQ,dQ1;
	int stack[L2];

	// escolho um spin, coloco na pilha e flipo ele
	i = FRANDOM*L2;
	stack[0] = i;
	sp = 1;
	oldspin = spin[i];
	//newspin = -spin[i];//(Qi + 1	+ (Q-2)*FRANDOM)%Q
	//newspin = spin[i] + (1 + (_Q-1)*FRANDOM)%_Q; // Q-1 ou Q-2 ?
	dQ1 = (_Q-1)*FRANDOM;
	dQ = (spin[i] + 1 + dQ1)%_Q; // testar se esta voltando pro Q original
	//newspin = spin[i] + dQ;
	newspin = dQ;
	if(newspin==_Q) {
		fprintf(stderr,"erro no newspin: %d %d %d (WOLFF)\n",newspin,dQ,dQ1);
		exit(-1);
	}
	spin[i] = newspin;
	//gnuplot_view(sp,spin);
	
	while(sp) {
		// pego um sitio da pilha
		// e diminuo o indice sp depois de pegar o valor
		current = stack[--sp];
		// confere os vizinhos
		for(j=0;j<4;j++) {
			nn = neigh[current][j];
			// posso juntar as duas condicoes. faz diferenca?
			if(spin[nn]==oldspin) {
	if(FRANDOM<exp_factor) {
		// coloco no proximo indice da pilha
		stack[sp++] = nn;
		// atualizo o estado do vizinho. dessa forma nao pode entrar
		// na pilha novamente ja que tem estado newspin
		spin[nn] = newspin;
		//gnuplot_view(sp,spin);
	}
			}
		}
	}
	//
	//gnuplot_view(sp,spin);
	//
	return;
}

/**************************************************************************
 * 			 Visualization	routine				  *
 *	compile with: -DGNUPLOT						  *
 *	use as: ./a.out | gnuplot					  *
 *************************************************************************/
void gnuplot_view(int tempo, int *s)
{
	int i,j,sitio;

	printf("set title \'tempo: %d \' \n",tempo);
	printf("set size square\n");
	printf("plot \'-\' matrix with image\n");

	for(j=0;j<L;j++)
	{
		for(i=0;i<L;i++)
		{
			sitio = i+j*L;
			printf(" %d",s[sitio]);
		}
		printf("\n"); 
	}

	printf("e\n pause 0.05\n"); 
	printf("\n\n"); 
	 
	return;
}
