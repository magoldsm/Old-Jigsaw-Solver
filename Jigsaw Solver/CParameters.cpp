#include "pch.h"
#include "CParameters.h"

using namespace std;

//#define LARGE_PUZZLE			// >= 100
//#define MEDIUM_PUZZLE			// >= 69
#define SMALL_PUZZLE			// < 69

constexpr auto APP = "Jigsaw Solver";

#ifdef SMALL_PUZZLE
#define	N_PARAMS	2

static int defaultP[] = { 940, 900 };
static int defaultM[] = { 2000, 2000 };
static int defaultMu[] = { 600, 600 };
static int defaultK3[] = { 707, 2000 };

#define	LAMBDA0		30
#define	LAMBDA1		10
#define	Q_1			900
#define	Q_2			200
#define	Q_2STAR		300
#define	Q_3			0
#endif // SMALL_PUZZLE

#ifdef MEDIUM_PUZZLE
#define	N_PARAMS	2

static int defaultP[] = { 940, 900 };
static int defaultM[] = { 3000, 2000 };
static int defaultMu[] = { 600, 600 };
static int defaultK3[] = { 1118, 2000 };

#define	LAMBDA0		20
#define	LAMBDA1		10
#define	Q_1			900
#define	Q_2			200
#define	Q_2STAR		300
#define	Q_3			0
#endif // MEDIUM_PUZZLE

#ifdef LARGE_PUZZLE
#define	N_PARAMS	1

static int defaultP[] = { 940 };
static int defaultM[] = { 2000 };
static int defaultMu[] = { 600 };
static int defaultK3[] = { 707 };

#define	LAMBDA0		20
#define	LAMBDA1		8
#define	Q_1			976
#define	Q_2			240
#define	Q_2STAR		300
#define	Q_3			29
#endif // LARGE_PUZZLE

CParameters::CParameters()
{
	sgOrder = ::GetProfileInt(APP, "sgOrder", 2);			// 5 for raw data
	sgWindow = ::GetProfileInt(APP, "sgWindow", 10);		// 50 for raw data
	lambda0 = ::GetProfileInt(APP, "lambda0", LAMBDA0);
	lambda1 = ::GetProfileInt(APP, "lambda1", LAMBDA1);
	alpha = ::GetProfileInt(APP, "alpha", 2);
	beta = ::GetProfileInt(APP, "beta", 5);
	gamma = ::GetProfileInt(APP, "gamma", 5);
	nu = ::GetProfileInt(APP, "nu", 3);
	rho = ::GetProfileInt(APP, "rho", 333) / 1000.0;
	C1 = ::GetProfileInt(APP, "C1", 1000);
	C2 = ::GetProfileInt(APP, "C2", 2);
	K1 = ::GetProfileInt(APP, "K1", 15);
	K2 = ::GetProfileInt(APP, "K2", 4);
	K4 = ::GetProfileInt(APP, "K3", 500)/1000.0;
	epsilon = ::GetProfileInt(APP, "epsilon", 100) / 1000000.0;
	eta1 = ::GetProfileInt(APP, "eta1", 1000) / 1000.0;
	eta2 = ::GetProfileInt(APP, "eta2", 1500) / 1000.0;
	Q1 = ::GetProfileInt(APP, "Q1", Q_1)/1000.0;
	Q2 = ::GetProfileInt(APP, "Q2", Q_2) / 1000.0;
	Q2Star = ::GetProfileInt(APP, "Q2Star", Q_2STAR) / 1000.0;
	Q3 = ::GetProfileInt(APP, "Q3", Q_3) / 1000.0;
	K2 = ::GetProfileInt(APP, "K2", 4);
	jMax = ::GetProfileInt(APP, "jMax", 3);

	int nParams = ::GetProfileInt(APP, "nParams", N_PARAMS);
	seq.resize(nParams);

	for (int i = 0; i < nParams; i++)
	{
		seq[i].p0 = ::GetProfileInt(APP, (string("p0_") + to_string(i)).data(), defaultP[i]) / 1000.0;
		seq[i].m0 = ::GetProfileInt(APP, (string("m0_") + to_string(i)).data(), defaultM[i]) / 1000.0;
		seq[i].mu0 = ::GetProfileInt(APP, (string("mu0_") + to_string(i)).data(), defaultMu[i]) / 1000.0;
		seq[i].K3 = ::GetProfileInt(APP, (string("K3_") + to_string(i)).data(), defaultK3[i]) / 1000.0;
	}
}


CParameters::~CParameters()
{
}
