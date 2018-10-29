#define _CRT_SECURE_NO_WARNINGS
#include "random.h"
#include <stdio.h>
#include <math.h>
#include <Windows.h>

#define N_Part 300
#define Temp 120
#define SIGMA_INITIAL 3.624
#define E_INITIAL 159.5197
#define R_CUT_A 12.5
double BOX_LENGHT_X = 25.0;
double Y_X = 1.0;
double Z_X = 1.0;


#define OUTPUT_FREQUENCY 1000000
#define N_STEPS 100000000

#define LONG_INT // LONG_INT LONG_LONG_INT FLOAT DOUBLE
#define CUBIC // CUBIC RECTANGULAR
int AcceptStep = 0;
double sigmaSquared;
double epsilonReduce;

#ifdef LONG_INT
double R_CUT = pow(2.0, 64) * R_CUT_A / BOX_LENGHT_X * R_CUT_A / BOX_LENGHT_X;
int coordx[N_Part];
int coordy[N_Part];
int coordz[N_Part];
int dR[3];
int dRMax = pow(2.0, 31) / 10.0;
double MAXR = pow(2, 32);
#endif // LONG_INT
#ifdef LONG_LONG_INT
double R_CUT = pow(2, 128) * R_CUT_A / BOX_LENGHT_X * R_CUT_A / BOX_LENGHT_X;
long long int coordx[N_Part];
long long int coordy[N_Part];
long long int coordz[N_Part];
long long int dR[3];
long long int dRMax = pow(2.0, 63) / 10.0;
double MAXR = pow(2, 64);
#endif // LONG_LONG_INT
#ifdef FLOAT
float check_pbc2(float r)
{
	r -= int(((r<0) ? -0.5 : 0.5) + r);
	return r;
}
float check_pbc(float r)
{
	r = abs(r);
	if (r > 0.5)
	{
		r -= 1.0;
	}
	return r;
}
double R_CUT = R_CUT_A / BOX_LENGHT_X * R_CUT_A / BOX_LENGHT_X;
float coordx[N_Part];
float coordy[N_Part];
float coordz[N_Part];
float dR[3];
float dRMax = 0.5 / 10.0;
double MAXR = 1.0;
#endif // FLOAT
#ifdef DOUBLE
double check_pbc2(double r)
{
	r -= int(((r<0) ? -0.5 : 0.5) + r);
	return r;
}
double check_pbc(double r)
{
	r = abs(r);
	if (r > 0.5)
	{
		r -= 1.0;
	}
	return r;
}
double R_CUT = R_CUT_A / BOX_LENGHT_X* R_CUT_A / BOX_LENGHT_X;
double coordx[N_Part];
double coordy[N_Part];
double coordz[N_Part];
double dR[3];
double dRMax = 0.5 / 10.0;
double MAXR = 1.0;
#endif // DOUBLE
void Init_coord()
{
	for (int i = 0; i < N_Part; i++)
	{
		coordx[i] = rand_half() * MAXR;
		coordy[i] = rand_half() * MAXR;
		coordz[i] = rand_half() * MAXR;
	}
	return;
}
#ifdef LONG_INT
double Eng(int N, int coordNx, int coordNy, int coordNz)
#endif
#ifdef LONG_LONG_INT
double Eng(int N, long long int coordNx, long long int coordNy, long long int coordNz)
#endif
#ifdef FLOAT
double Eng(int N, float coordNx, float coordNy, float coordNz)
#endif
#ifdef DOUBLE
double Eng(int N, double coordNx, double coordNy, double coordNz)
#endif
{
	double R;
	double eng = 0;
	double buf;
	for (int i = 0; i < N_Part; i++)
	{
		if (N != i)
		{
#ifdef CUBIC
#if defined(LONG_INT) || defined(LONG_LONG_INT)
			buf = coordx[i] - coordNx;
			R = buf * buf;
			buf = coordy[i] - coordNy;
			R += buf * buf;
			buf = coordz[i] - coordNz;
			R += buf * buf;
#else 
			buf = check_pbc(coordx[i] - coordNx);
			R = buf * buf;
			buf = check_pbc(coordy[i] - coordNy);
			R += buf * buf;
			buf = check_pbc(coordz[i] - coordNz);
			R += buf * buf;
#endif // Coordinate storage type
#else 
#endif // CUBIC
#ifdef RECTANGULAR
#if defined(LONG_INT) || defined(LONG_LONG_INT)
			buf = coordx[i] - coordNx;
			R = buf * buf;
			buf = (coordy[i] - coordNy) * Y_X;
			R += buf * buf;
			buf = (coordz[i] - coordNz) * Z_X;
			R += buf * buf;
#else 
			buf = check_pbc(coordx[i] - coordNx);
			R = buf * buf;
			buf = check_pbc(coordy[i] - coordNy) * Y_X;
			R += buf * buf;
			buf = check_pbc(coordz[i] - coordNz) * Z_X;
			R += buf * buf;
#endif // Coordinate storage type
#else 
#endif // RECTANGULAR
			if (R < R_CUT)
			{
				R = sigmaSquared / R;
				R = R * R * R;
				eng += R * R - R;
			}
		}
	}
	return eng * epsilonReduce;
}

double Energy()
{
	double R;
	double eng = 0;
	double buf;
	for (int i = 0; i < N_Part; i++)
	{
		for (int j = 0; j < i; j++)
		{
#ifdef CUBIC
#if defined(LONG_INT) || defined(LONG_LONG_INT)
			buf = coordx[i] - coordx[j];
			R = buf * buf;
			buf = coordy[i] - coordy[j];
			R += buf * buf;
			buf = coordz[i] - coordz[j];
			R += buf * buf;
#else 
			buf = check_pbc(coordx[i] - coordx[j]);
			R = buf * buf;
			buf = check_pbc(coordy[i] - coordy[j]);
			R += buf * buf;
			buf = check_pbc(coordz[i] - coordz[j]);
			R += buf * buf;
#endif // Coordinate storage type
#else 
#endif // CUBIC
#ifdef RECTANGULAR
#if defined(LONG_INT) || defined(LONG_LONG_INT)
			buf = coordx[i] - coordx[j];
			R = buf * buf;
			buf = (coordy[i] - coordy[j]) * Y_X;
			R += buf * buf;
			buf = (coordz[i] - coordz[j]) * Z_X;
			R += buf * buf;
#else 
			buf = check_pbc(coordx[i] - coordx[j]);
			R = buf * buf;
			buf = check_pbc(coordy[i] - coordy[j]) * Y_X;
			R += buf * buf;
			buf = check_pbc(coordz[i] - coordz[j]) * Z_X;
			R += buf * buf;
#endif // Coordinate storage type
#else 
#endif // RECTANGULAR
			if (R < R_CUT)
			{
				R = sigmaSquared / R;
				R = R * R * R;
				eng += R * R - R;
			}
		}
	}
	return eng * epsilonReduce;
}
void random_dR()
{
	dR[0] = rand_half() * dRMax;
	dR[1] = rand_half() * dRMax / Y_X;
	dR[2] = rand_half() * dRMax / Z_X;
	return;
}
void stepUpdate()
{
	if ((double)AcceptStep / OUTPUT_FREQUENCY > 0.5)
	{
		dRMax *= 1.05;
	}
	else
	{
		dRMax /= 1.05;
	}
	AcceptStep = 0;
}


int main()
{
	
#ifdef LONG_INT
	sigmaSquared = SIGMA_INITIAL / BOX_LENGHT_X * pow(2, 32)*SIGMA_INITIAL / BOX_LENGHT_X * pow(2, 32);
#endif
#ifdef LONG_LONG_INT
	sigmaSquared = SIGMA_INITIAL / BOX_LENGHT_X * pow(2, 64)*SIGMA_INITIAL / BOX_LENGHT_X * pow(2, 64);
#endif
#ifdef FLOAT
	sigmaSquared = SIGMA_INITIAL / BOX_LENGHT_X *SIGMA_INITIAL / BOX_LENGHT_X;
#endif
#ifdef DOUBLE
	sigmaSquared = SIGMA_INITIAL / BOX_LENGHT_X *SIGMA_INITIAL / BOX_LENGHT_X;
#endif
	epsilonReduce = 4.0 * E_INITIAL / Temp;

	double eng = 0;
	double engOld;
	double engNew;
	double engAverage = 0;
	int chosenN;

	Init_coord();
	double tail = 3.14159265358979323846 * 2.0 / 3.0 / pow(R_CUT_A*2.0, 3)*epsilonReduce*pow(SIGMA_INITIAL, 3)*(1 / 3 * pow(R_CUT_A / SIGMA_INITIAL, -9) - pow(R_CUT_A / SIGMA_INITIAL, -3)) * N_Part * N_Part;

	eng = Energy() + tail;
	int start_time = clock();

	for (int nSteps = 0; nSteps < N_STEPS; nSteps++)
	{
		chosenN = rand_zero_one() *N_Part;
		random_dR();
		// for different interaction can be converted to cycle
		engOld = Eng(chosenN, coordx[chosenN], coordy[chosenN], coordz[chosenN]);
#if defined(LONG_INT) || defined(LONG_LONG_INT)
		engNew = Eng(chosenN, coordx[chosenN] + dR[0], coordy[chosenN] + dR[1], coordz[chosenN] + dR[2]);
#else
		engNew = Eng(chosenN, check_pbc2(coordx[chosenN] + dR[0]), check_pbc2(coordy[chosenN] + dR[1]), check_pbc2(coordz[chosenN] + dR[2]));
#endif
		if (exp(engOld - engNew) > rand_zero_one())
		{
			eng = eng + engNew - engOld;
#if defined(LONG_INT) || defined(LONG_LONG_INT)
			coordx[chosenN] += dR[0];
			coordy[chosenN] += dR[1];
			coordz[chosenN] += dR[2];
#else
			coordx[chosenN] = check_pbc2(coordx[chosenN] + dR[0]);
			coordy[chosenN] = check_pbc2(coordy[chosenN] + dR[1]);
			coordz[chosenN] = check_pbc2(coordz[chosenN] + dR[2]);
#endif
			AcceptStep++;
		}
		engAverage += eng;
		if (nSteps == OUTPUT_FREQUENCY)
		{
			eng = Energy() + tail;
		}
		if (nSteps % OUTPUT_FREQUENCY == 0 && nSteps != 0)
		{
			stepUpdate();
			printf("%d\t%d\t%lf\n", nSteps, clock() - start_time, engAverage / OUTPUT_FREQUENCY);
			engAverage = 0;
		}

	}
	system("pause");
	return 0;
}

