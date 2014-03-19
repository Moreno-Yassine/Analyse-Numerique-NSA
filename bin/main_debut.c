#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "von_neumann.h"
#include "aes.h"
#include "mersenne_twister.h"

#define ARRAY_MAX_SIZE 1000
#define OLDRAND_MAX 2147483647
#define TAILLE 1024

static int next;

int rdtsc()
{
	// cette fonction suivante cree un warning : c'est normal.
	__asm__ __volatile__("rdtsc");
}

void oldinit_rand(int seed)
{
	next = seed;
}

int oldrand()
{
	next = next * 1103515245 + 12345;
	return (unsigned int)(next % OLDRAND_MAX);
}

// FONCTION doubleFrequency
double doubleFrequency(int n, unsigned int* resultat, int size)
{
	double somme = 0;
		for (int k=0;k<size;k++)
		{
			int mask = 1;
			for (int bit = 0; bit<n ;bit++)
			{
				if((resultat[k]) & (mask))
				{
					somme += 1;
				}
				else
				{
					somme -= 1;
				}
				mask = mask* 2 ;
			}
		}		
	double sommeobs =  abs(somme)/sqrt(size*n);	
	return sommeobs;
}

// FONCTION doubleRuns
double doubleRuns(int n, unsigned int* resultat, int size)
{
	double somme = 0;
	double pi = 0;
		for (int k=0;k<size;k++)
		{
			int mask = 1;
			for (int bit = 0; bit<n ;bit++)
			{
				if((resultat[k]) & (mask))
				{
					somme += 1;
				}
				mask = mask* 2 ;
			}
		}		
	pi = somme/(size*n);	
	// PRE TEST
	double tau = 2/sqrt(size*n);
	
	if (abs(pi-0.5)>=tau)
	{
	return 0;
	}
	else
	{
		double vobs = 1;
		int previous = 0;
		int current = 0;
		for (int k=0;k<size;k++)
		{
			int mask = 1;
			for (int bit = 0; bit<n ;bit++)
			{
				if((resultat[k]) & (mask))
				{
					current = 1;
				}
				else
				{
					current = 0;
				}
				if (current!=previous)
				{
					vobs++;
				}
				previous = current;
				mask = mask* 2 ;
			}
		}	
		double num = abs(vobs - 2*size*n*pi*(1-pi));
		double denum = 2*sqrt(2*size*n)*pi*(1-pi);
		double pvaleur = erfc (num/denum);
		return pvaleur;
	}
	
	

}

int main()
{
	word16 x=1111; // nombre entre 1000 et 9999 pour Von Neumann
	struct mt19937p mt; // Pour Mersenne-Twister
	int tmp = rand(), seed; // Pour Mersenne-Twister
	u32 Kx[NK], Kex[NB*NR], Px[NB]; // pour l'AES

	int output_rand; // sortie du rand du C	
	word32 output_AES; // sortie pour l'AES
	word16 output_VN; // sortie pour pour Von Neumann
	word32 output_MT; // sortie pour Mersenne-Twister

               
	// initialisation des graines des generateurs

	srand(rdtsc());  // rand du C 
	seed = rand();
	oldinit_rand(seed);
	sgenrand(time(NULL)+(tmp), &mt); // Mersenne-Twister
	// Initialisation de la clé et du plaintext pour l'AES 
	// 45 est un paramètre qui doit changer à chaque initialisation
	init_rand(Kx, Px, NK, NB, 45);
	KeyExpansion(Kex,Kx); // AES : sous-clefs


	// sorties des generateurs	
	output_rand = oldrand(); // rand du C
	output_VN = Von_Neumann(&x); // Von Neumann
	output_MT = genrand(&mt); // Mersenne-Twister
	output_AES = AES(Px, Kex); // AES
	
	FILE * pfile;
	pfile= fopen ("SortieTest.txt","w");
	for (int cycle =1 ;cycle<21;cycle++)
	{
		fprintf(pfile,"- Cycle %i - \n",cycle);
		unsigned int pFaible [TAILLE];
		unsigned int pFort [TAILLE];
		unsigned int vn [TAILLE];
		unsigned int mtwister [TAILLE];
		unsigned int aes [TAILLE];
		// Generation de séquence
		for (int bloc=0;bloc<TAILLE;bloc ++)
		{
			 pFaible [bloc] = oldrand()&0x0f;
			 pFort [bloc]= (oldrand()>>27)&0x0f;
			 vn [bloc] = Von_Neumann(&x);
			 mtwister [bloc] = genrand(&mt);
			 aes [bloc] = AES(Px, Kex);
		}
		// TESTS MONOBIT
		double outFaible = 0;
		double outFort = 0;
		double outvn = 0;
		double outmt = 0;
		double outAES = 0;

			outFaible = doubleFrequency (4,pFaible,TAILLE);
			outFort = doubleFrequency (4,pFort,TAILLE);
			outvn = doubleFrequency (13,vn,TAILLE); // D'après l'histogramme de VN
			outmt = doubleFrequency (32,mtwister,TAILLE); // d'après l'histogramme de MT
			outAES = doubleFrequency (32,aes,TAILLE);

		fprintf(pfile,"- TESTS MONOBIT - \n");	
		fprintf(pfile,"- Test Monobit rand Faible : %f \n",erfc(outFaible/sqrt(2)));
		fprintf(pfile,"- Test Monobit rand Fort : %f \n",erfc(outFort/sqrt(2)));
		fprintf(pfile,"- Test Monobit Von Neumann : %f \n",erfc(outvn/sqrt(2)));
		fprintf(pfile,"- Test Monobit Mersenne - Twister : %f \n",erfc(outmt/sqrt(2)));
		fprintf(pfile,"- Test Monobit AES : %f \n",erfc(outAES/sqrt(2)));	
		fprintf(pfile, "\n");

			outFaible = doubleRuns (4,pFaible,TAILLE);
			outFort = doubleRuns (4,pFort,TAILLE);
			outvn = doubleRuns (13,vn,TAILLE); // D'après l'histogramme de VN
			outmt = doubleRuns (32,mtwister,TAILLE); // d'après l'histogramme de MT
			outAES = doubleRuns (32,aes,TAILLE);
			
		fprintf(pfile,"- TESTS RUNS - \n");	
		fprintf(pfile,"- Test RUNS rand Faible : %f \n",outFaible);
		fprintf(pfile,"- Test RUNS rand Fort : %f \n",outFort);
		fprintf(pfile,"- Test RUNS Von Neumann : %f \n",outvn);
		fprintf(pfile,"- Test RUNS Mersenne - Twister : %f \n",outmt);
		fprintf(pfile,"- Test RUNS AES : %f \n",outAES);
		
		fprintf(pfile, "\n");
		fprintf(pfile, "\n");
			
	}
	
	
	
	
	// affichage
	printf("- Generation de nombres aleatoires -\n");
	printf("rand du C : %u \n",output_rand); 
	printf("Von Neumann : %u\n",output_VN);
	printf("Mersenne Twister : %u\n",output_MT);
	printf("AES : %u\n",output_AES);


	return 0;
}
