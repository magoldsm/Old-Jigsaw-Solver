#pragma once
class CParameters
{
public:
	CParameters();
	~CParameters();

	int sgOrder;
	int sgWindow;
	int lambda0;
	int lambda1;
	int alpha;
	int beta;
	int gamma;
	int nu;
	double rho;
	int jMax;
	int C1;
	int C2;
	double K1;
	double K2;
	double K4;
	double epsilon;

	struct ParameterSequence
	{
		double	p0;
		double	m0;
		double	mu0;
		double	K3;
	};

	std::vector<ParameterSequence>	seq;

};

static CParameters params;				// The one and only!
