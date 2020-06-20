#include "pch.h"
#include "CPiece.h"
#include "SignatureSimilarity.h"
#include "CParameters.h"
#include "Utilities.h"
#include "Jigsaw Solver.h"
#include "Lock.h"

using namespace Eigen;
using namespace std;

using VectorXb = Matrix<bool, -1, 1>;

/*
function new_points = transf(points, trans, cm)

% This function applies a rigid motion to inputed points
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'points':   This should be n points represented as an n-by-2 matrix,
			each row of which specifies a point in R^2.

'trans':    This input should be a matrix of the form [theta a b]
			where these parameters specify a rigid motion consisting
			of a translation by [a b] and a rotation of theta radians
			around the point cm.

'cm':       The point in R^2, represented by a 1-by-2 matrix, about
			which the rotation will take place.

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'new_points':   This output gives the transformed points as an n-by-2
			matrix, each row of which specifies a point in R^2.

%--------------------------------------------------------------------
%}
*/

inline static Matrix<double, 1, 2>
TransformAboutCM(const Vector2d point, const GTransform& trans, Vector2d cm)
{
	Vector2d result;

	double c = cos(trans.theta);
	double s = sin(trans.theta);

	result[0] = (point(0) * c - point(1) * s) + trans.dx + cm[0];
	result[1] = (point(0) * s + point(1) * c) + trans.dy + cm[1];

	return result;
}

/*
% This function performs the locking algorithm described in [1]. It
% is intended for use by the Assemble() function.
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'g_0':      This input should be a matrix of the form [theta a b]
			where these parameters specify a rigid motion as
			described in [2] (a translation by [a b] and a rotation
			by theta radians around the origin).

'C_Delta':  This should be a discretized planar curve of n points
			represented as an n-by-2 matrix, each row of which
			specifies a point in R^2. The curve is assumed to be
			closed; no repetition of points is necessary.

'Ctilde_Delta': This should be a discretized planar curve of n points
			represented as an n-by-2 matrix, each row of which
			specifies a point in R^2. The curve is assumed to be
			closed; no repetition of points is necessary.

'K_1', 'K_2', 'K_3', 'K_4', 'epsilon', 'nu', 'rho', 'j_max':   These inputs
			should be doubles that give the parameters as described
			in [1].

'plotter':  In this bool, value of 1 results in the plotting of a
			visualization of the piece locking whereas a value of 0
			does not.

'fh':       If this input is supplied and plotter = 1, the
			visualization is plotted on the figure with handle 'fh'.

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'g_lock'	: This is the resulting transformation as described in
			[1].

'D_Delta_3_Ics', 'Dtilde_Delta_3_Ics', 'D_Delta_2_Ics',
'Dtilde_Delta_2_Ics': These columns contain the
			indices (within C_Delta and Ctilde_Delta) of
			the points in the sets as defined in [1].

'xph':      This output is assigned only if 'fh' is inputted and
			plotter = 1, in which case it contains the handle of an
			object to be later deleted by the Assemble() function.

%--------------------------------------------------------------------
%}
*/

void
Lock(const GTransform& g0, const Curve& CDelta, const Curve& CtildeDelta, GTransform& gLock, double K3, Indices& D_Delta_3_Ics, Indices& Dtilde_Delta_3_Ics, Indices& D_Delta_2_Ics, Indices& Dtilde_Delta_2_Ics)
{
	int nC = (int)CDelta.rows();
	int nCtilde = (int)CtildeDelta.rows();
	Curve CDeltam1;
	vector<Curve> EtildeDeltaK1;
	Curve EDeltaK1;
	double thetaj = 0.0;
	Vector2d cj(0.0, 0.0);

	circShift(CDelta, CDeltam1, -1);

	// Equation 5.3

	double dStar = ((CDelta.col(0) - CDeltam1.col(0)).array().square() + (CDelta.col(1) - CDeltam1.col(1)).array().square()).sqrt().sum() / nC;
	double dStarK1 = params.K1*dStar;

	// Equations 5.4
	//
	// bNearby is a boolean matrix.  It has CtildeDelta.rows() rows by CDelta.rows() columns.  Thus, each row corresponds
	// to an element (ztilde) in CtildeDelta) and each column corresponds to an element in CDelta.
	//
	// An entry is true if || ztilde - g . z || < dStar * K1.  That is if the points are "close enough" to interact in the
	// locking alogorithm
	// 

	Matrix<bool, -1, -1> bNearby;
	bNearby.resize(CtildeDelta.rows(), CDelta.rows());

	FOR_START(c1, 0, nC)
		bNearby.col(c1) = (CtildeDelta.rowwise() - TransformAboutCM(CDelta.row(c1), g0, Vector2d(0, 0))).rowwise().norm().array() < dStarK1;
	FOR_END

		// The goal is now to create 2 sets: EDeltaK1 and EtildeDeltaK1.  EtildeDeltaK1 is the set of points
		// Tots is the number of trues in each column of bNearby.   nEDeltas
		Matrix<Index, 1, -1> tots = bNearby.colwise().count();
	int nEDeltas = (int)(tots.array() > 0).count();

	if (nEDeltas == 0)
	{
		gLock = GTransform();
		return;
	}

	EtildeDeltaK1.resize(nEDeltas);
	EDeltaK1.resize(nEDeltas, 2);
	int iEtildeDeltaK1 = 0;

	for (int c1 = 0; c1 < nC; c1++)
	{
		if (tots[c1] != 0)
		{
			EDeltaK1.row(iEtildeDeltaK1) = CDelta.row(c1);
			Curve& c = EtildeDeltaK1[iEtildeDeltaK1++];
			c.resize(tots[c1], 2);
			int idx = 0;
			for (int r1 = 0; r1 < nCtilde; r1++)
			{
				if (bNearby(r1, c1))
				{
					c.row(idx++) = CtildeDelta.row(r1);
				}
			}
		}
	}

	// Calculate perturbation constants.  Equations 5.5 and 5.6

	RowVector2d zCM = EDeltaK1.colwise().sum() / nEDeltas;
	RowVectorXd normSquared = (EDeltaK1.rowwise() - zCM).rowwise().squaredNorm();
	double r2 = normSquared.sum();
	double rINF = normSquared.array().sqrt().maxCoeff();

	// Convert g_0 from the form used in Assemble() (rotation around the
	// origin) to the form used in Lock() (rotation around the center of
	// mass

	double c = cos(g0.theta);
	double s = sin(g0.theta);
	Matrix2d rot;
	rot << c, -s, s, c;

	Vector2d zCMt = zCM.transpose();

	MatrixXd translation = rot * zCMt - zCMt + Vector2d(g0.dx, g0.dy);

	GTransform gj(g0.theta, translation(0), translation(1));
	RowVector2d wj = zCMt + translation;

	// Perform iterative pertubation

	double dj = dStarK1;
	double ThetaJm1 = 0.0;
	Vector2d cJm1(0.0, 0.0);
	double dStarK4 = params.K4 * dStar;

	for (int j = 0; j < params.jMax; j++)
	{
		double tauTotj = 0.0;
		RowVector2d fTotj(0.0, 0.0);
		VectorXd Aj(nEDeltas);

		for (int c2 = 0; c2 < nEDeltas; c2++)
		{
			RowVector2d EDeltaK1Transformed = TransformAboutCM(EDeltaK1.row(c2), gj, zCM);
			MatrixXd diff = EtildeDeltaK1[c2].rowwise() - EDeltaK1Transformed;
			VectorXd diffNorm = diff.rowwise().norm();
			Aj[c2] = diffNorm.minCoeff();

			Vector2d fj;

			if (Aj[c2] >= dStarK4)
			{
				Vector2d fj = (diff.array().colwise() / (diffNorm.array().pow(params.nu + 1) + params.epsilon * diffNorm.array())).colwise().sum();
				fTotj.array() += fj.array();
				RowVector2d t = EDeltaK1Transformed - wj;
				tauTotj += t[0] * fj[1] - t[1] * fj[0];		// 
			}
		}

		// Equation 5.10

		double dAvj = Aj.mean();
		qsort(Aj.data(), Aj.rows(), sizeof(double), [](const void* p1, const void* p2) -> int
		{
			return (int)(*(double*)p1 - *(double*)p2);
		});

		double dMedj = Aj.rows() & 1 ? Aj[Aj.rows() / 2] : (Aj[Aj.rows() / 2 - 1] + Aj[Aj.rows() / 2]) / 2;

		// Terminate the algorithm if fit is poor and getting worse.

		if (dAvj > dj && dMedj >= params.seq[j].K3*dStar)
		{
			gj.theta -= thetaj;
			gj.dx -= cj[0];
			gj.dy -= cj[1];
			// ***BUGBUG*** Plotting code
			break;
		}

		// Calculate new transformation.  Equation 5.12

		double deltaj = params.rho * dMedj / max(fTotj.norm() / nEDeltas, rINF*fabs(tauTotj) / r2);

		// Equations 5.11

		thetaj = deltaj * tauTotj / r2;
		cj = deltaj * fTotj / nEDeltas;

		// Update key quantities - Step 7

		gj.theta += thetaj;
		gj.dx += cj[0];
		gj.dy += cj[1];
		wj += cj;
		dj = dAvj;

		// ***BUGBUG*** Plotting code

		// Step 7, termination condition

		if ((thetaj*ThetaJm1 < 0 && cj[0] * cJm1[0] < 0 && cj[1] * cJm1[1] < 0))
			break;

		ThetaJm1 = thetaj;
		cJm1 = cj;
	}

	// Find indices of relevant sets

	double K2dStar = params.K2 * dStar;
	double K3dStar = K3 * dStar;
	//Matrix<bool, -1, 1> D2Sel(nCtilde);
	//Matrix<bool, -1, 1> D3Sel(nCtilde);

	//D2Sel.setConstant(false);
	//D3Sel.setConstant(false);

	for (int c2 = 0; c2 < nC; c2++)
	{
		VectorXd diff = (CtildeDelta.rowwise() - TransformAboutCM(CDelta.row(c2), gj, zCM)).rowwise().norm();
		VectorXb test = diff.array() < K2dStar;

		if (test.any())
		{
			for (int i = 0; i < nCtilde; i++)
				if (test[i])
					Dtilde_Delta_2_Ics.push_back(i);
			D_Delta_2_Ics.push_back(c2);
		}

		test = diff.array() < K3dStar;
		if (test.any())
		{
			for (int i = 0; i < nCtilde; i++)
				if (test[i])
					Dtilde_Delta_3_Ics.push_back(i);
			D_Delta_3_Ics.push_back(c2);
		}
	}

	auto icmp = [](const void* p1, const void* p2) { return *(int*)p1 - *(int*)p2; };
	qsort(Dtilde_Delta_2_Ics.data(), Dtilde_Delta_2_Ics.size(), sizeof(int), icmp);
	qsort(Dtilde_Delta_3_Ics.data(), Dtilde_Delta_3_Ics.size(), sizeof(int), icmp);

	auto ucmp = [](int i, int j) { return i == j; };
	unique(Dtilde_Delta_2_Ics.begin(), Dtilde_Delta_2_Ics.end(), ucmp);
	unique(Dtilde_Delta_3_Ics.begin(), Dtilde_Delta_3_Ics.end(), ucmp);

	// ***BUGBUG*** Plotting code

	// Convert gj to the form used in SignatureSimilarity()

	//g_lock = [g_j(1), (g_j(2:3)' + z_cm' - [cos(g_j(1)), -sin(g_j(1)); sin(g_j(1)) cos(g_j(1))] * z_cm')'];

	c = cos(gj.theta);
	s = sin(gj.theta);
	rot << c, -s, s, c;

	translation = Vector2d(gj.dx, gj.dy) + zCMt - rot * zCMt;

	gLock = GTransform(gj.theta, translation(0), translation(1));
	int yyyy = 0;
}
