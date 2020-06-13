#include "pch.h"
#include "CPiece.h"
#include "SignatureSimilarity.h"
#include "CParameters.h"
#include "Utilities.h"
#include "Jigsaw Solver.h"
#include "Lock.h"

using namespace Eigen;
using namespace std;
using namespace cv;

#pragma once



CPScore::CPScore(size_t sz) : m_Size(sz)
{
	m_ArcScores = new MatrixXd[sz*sz];
}


CPScore::~CPScore()
{
	if (m_ArcScores) delete[] m_ArcScores;
}

void CPScore::SetSize(size_t sz)
{
	if (m_ArcScores) delete[] m_ArcScores;
	m_Size = sz;
	m_ArcScores = new MatrixXd[sz*sz];
}

void CPScore::Display()
{
	cout << "    ";
	for (size_t i = 0; i < m_Size; i++)
	{
		cout << setw(5) << i << " ";
	}
	cout << endl;
	for (size_t i = 0; i < m_Size; i++)
	{
		cout << setw(2) << i << "  ";
		for (size_t j = 0; j < m_Size; j++)
		{
			MatrixXd& arcscore = (*this)(i, j);
			if (arcscore.cols() == 0 && arcscore.rows() == 0)
				cout << "   [] ";
			else
			{
				char buff[100];
				sprintf_s(buff, 100, "%dx%d ", (int) arcscore.rows(), (int) arcscore.cols());
				cout << setw(6) << buff;
			}
		}
		cout << endl;
	}
}

void CPScore::Display(size_t nRow, size_t nCol)
{
	MatrixXd& arcscore = (*this)(nRow, nCol);
	for (int i = 0; i < arcscore.rows(); i++)
	{
		for (int j = 0; j < arcscore.cols(); j++)
		{
			cout << setprecision(4) << setw(8) << arcscore(i, j);
		}
		cout << endl;
	}
}


// Equation numbers are from paper "Extensions of Invariant Signatures for Object Recognition

static double
_SignatureSimilarity(const Curve& sig1, const Curve& sig2, double D)
{
	int sz1 = (int) sig1.size();
	int sz2 = (int) sig2.size();

	VectorXd h = VectorXd::Zero(sz1);

	// Equations 3.4, 3.5

	for (int i = 0; i < sz1; i++)
	{
		Point2d t = sig1[i];
		VectorXd d = ((sig2.kappa.array() - t.x).array().square() + (sig2.kappas.array()-t.y).array().square()).sqrt();
		for (int j = 0; j < sz2; j++)
		{
			if (d[j] < D)
			{
				h[i] += 1.0 / (pow((d[j] / (D - d[j])), params.gamma) + params.epsilon);
			}
		}
	}

	// Equation 3.8

	VectorXd pPoint = h.array() / (h.array() + params.C1);

	// Equation 3.10

	VectorXd kPow = sig1.kappa.array().pow(params.alpha);
	VectorXd Prod = (pPoint.array() * kPow.array());
	double num = Prod.sum();
	double den = kPow.sum();
	return num/den;

}


double
SignatureSimilarity(const Curve& sig1, const Curve& sig2, double D)
{
	double pHat1 = _SignatureSimilarity(sig1, sig2, D);
	double pHat2 = _SignatureSimilarity(sig2, sig1, D);

	return min(pHat1, pHat2);
}

/*
% This function applies a rigid motion to inputted points
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'points':   This should be n points represented as an n-by-2 matrix,
			each row of which specifies a point in R^2.

'trans':    This input should be a matrix of the form [theta a b]
			where these parameters specify a rigid motion as
			described in [2] (a translation by [a b] and a rotation
			by theta radians around the origin).

%--------------------------------------------------------------------

%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'result':   This output gives the transformed points as an n-by-2
			matrix, each row of which specifies a point in R^2.

%--------------------------------------------------------------------
*/

Curve
TransformCurve(const Curve& points, const GTransform& trans)
{
	Curve result;

	double c = cos(trans.theta);
	double s = sin(trans.theta);

	result.x = (points.x*c - points.y*s).array() + trans.dx;
	result.y = (points.x*s + points.y*c).array() + trans.dy;

	return result;
}

/*
	%--------------------------------------------------------------------
	 Compute the rigid motion angle, as described in
	 Sect. 3.3.of[2]
	%--------------------------------------------------------------------
	INPUTS
	%--------------------------------------------------------------------
	
	'curve1':   This should be a discretized planar curve of n points
				represented as an n - by - 2 matrix, each row of which
				specifies a point in R^2.
	
	'curve2' : This should be a discretized planar curve of m points
				represented as an m - by - 2 matrix, each row of which
				specifies a point in R^2.
	
	'signature1' : This should be a discrete approximation to the
				Euclidean signature of 'curve1', represented as an
				n - by - 2 matrix, each row of which specifies a point in
				R^2.
	
	'signature2' : This should be a discrete approximation to the
				Euclidean signature of 'curve2', represented as an
				m - by - 2 matrix, each row of which specifies a point in
				R^2.
	
	'beta' :	This real number specifies the power by which kappa_s
				weights an average as described in Sect. 3.3.
	
	%--------------------------------------------------------------------
	
	
	
	%--------------------------------------------------------------------
	OUTPUTS
	%--------------------------------------------------------------------
	
	'theta':    The resulting rigid motion angle, as described in
				Sect. 3.3.of[2].
	
	%--------------------------------------------------------------------
	%
*/

static double
Rigid_Motion_Angle(const Curve& curve1, const Curve& curve2, const Curve& signature1, const Curve& signature2, double beta)
{
	VectorXd dx1, dx2, dy1, dy2;
	dx1.resizeLike(curve1.x);
	dx2.resizeLike(curve2.x);
	dy1.resizeLike(curve1.y);
	dy2.resizeLike(curve2.y);

	int sz1 = (int)dx1.size();
	int sz2 = (int)dx2.size();

	dx1.block(1, 0, sz1 - 2, 1) = curve1.x.topRows(sz1 - 2) - curve1.x.bottomRows(sz1 - 2);
	dx2.block(1, 0, sz2 - 2, 1) = curve2.x.topRows(sz2 - 2) - curve2.x.bottomRows(sz2 - 2);
	dy1.block(1, 0, sz1 - 2, 1) = curve1.y.topRows(sz1 - 2) - curve1.y.bottomRows(sz1 - 2);
	dy2.block(1, 0, sz2 - 2, 1) = curve2.y.topRows(sz2 - 2) - curve2.y.bottomRows(sz2 - 2);

	dx1[0] = dx1[1]; dx1[sz1 - 1] = dx1[sz1 - 2];
	dx2[0] = dx2[1]; dx2[sz2 - 1] = dx2[sz2 - 2];
	dy1[0] = dy1[1]; dy1[sz1 - 1] = dy1[sz1 - 2];
	dy2[0] = dy2[1]; dy2[sz2 - 1] = dy2[sz2 - 2];

	dx1.array() *= signature1.kappas.array();
	dy1.array() *= signature1.kappas.array();

	dx2.array() *= signature2.kappas.array();
	dy2.array() *= signature2.kappas.array();

	VectorXd ks1pow = signature1.kappas.array().pow(beta);
	VectorXd ks2pow = signature2.kappas.array().pow(beta);
	double ks1powsum = ks1pow.sum();
	double ks2powsum = ks2pow.sum();

	// The original code divides the X and Y componenets by ksXpowsum.  I believe this
	// is unnecessary and I've removed it.

	double angle1 = atan2((dy1.array() * ks1pow.array()).sum()/* / ks1powsum*/, (dx1.array() * ks1pow.array()).sum()/* / ks1powsum*/);
	double angle2 = atan2((dy2.array() * ks2pow.array()).sum()/* / ks2powsum*/, (dx2.array() * ks2pow.array()).sum()/* / ks2powsum*/);

	return fmod(angle2 - angle1 + 2 * PI, 2 * PI);
}

/*
% This function reconstructs a rigid motion translation as described
% in Sect. 3.3 of [2].
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'curve1':   This should be a discretized planar curve of n points
			represented as an n-by-2 matrix, each row of which
			specifies a point in R^2.

'curve2':   This should be a discretized planar curve of m points
			represented as an m-by-2 matrix, each row of which
			specifies a point in R^2.

'signature1':   This should be a discrete approximation to the
			Euclidean signature of 'curve1', represented as an
			n-by-2 matrix, each row of which specifies a point in
			R^2.

'signature2':   This should be a discrete approximation to the
			Euclidean signature of 'curve2', represented as an
			m-by-2 matrix, each row of which specifies a point in
			R^2.

'theta':    This should be the angle of the rigid motion as
			outputted by Rigid_Motion_Angle().

'beta':     This real number specifies the power by which kappa_s
			weights an average as described in Sect. 3.3.

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'G.dx':     The resulting rigid motion horizontal translation as
			described in Sect. 3.3 of [2].

'G.dy':     The resulting rigid motion vertical translation as
			described in Sect. 3.3 of [2].

%--------------------------------------------------------------------
%}
*/

static void
Rigid_Motion_Translation(GTransform& G, const Curve& curve1, const Curve& curve2, const Curve& signature1, const Curve& signature2, double theta, double beta)
{
	VectorXd ks1pow = signature1.kappas.array().pow(beta);
	VectorXd ks2pow = signature2.kappas.array().pow(beta);
	double ks1powsum = ks1pow.sum();
	double ks2powsum = ks2pow.sum();

	// Calculate center of mass weighted by kappa_s

	double cm1x = (curve1.x.array() * ks1pow.array()).sum() / ks1powsum;
	double cm1y = (curve1.y.array() * ks1pow.array()).sum() / ks1powsum;
	double cm2x = (curve2.x.array() * ks2pow.array()).sum() / ks2powsum;
	double cm2y = (curve2.y.array() * ks2pow.array()).sum() / ks2powsum;

	G.dx = cm2x - cm1x * cos(theta) + cm1y * sin(theta);
	G.dy = cm2y - cm1x * sin(theta) - cm1y * cos(theta);
}

/*
% This function calculates the mu score of a collection of rigid
% motions as described in [2].
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'trans':    This should be a n-by-3 matrix each row of which
			represents a transformation as [theta a b] where these
			parameters specify a rigid motion as described in [2].

'D_x':      The characteristic horizontal distance used in the mu
			score calculation as described in [2].

'D_y':      The characteristic vertical distance used in the mu score
			calculation as described in [2].

'C_2':      The parameter used in the mu score calculation as
			described in [2].

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'mu':       This output gives resulting mu score as described in [2].

%--------------------------------------------------------------------
%}

% Calculate differences between angles
difang = zeros(size(trans, 1), size(trans, 1));
for c1 = 1:size(trans, 1)
	for c2 = c1 + 1:size(trans, 1)
		difang(c1, c2) = mod(abs(trans(c1, 1) - trans(c2, 1)), 2*pi);
		if(difang(c1, c2) > pi)
			difang(c1, c2) = 2*pi - difang(c1, c2);
		end
	end
end

% Calculate mu score
mu = 1-C_2*sum([max(max(difang))/pi  range(trans(:, 2:3))./[D_x  D_y]], 2);
if mu < 0
	mu = 0;
end

*/

static GTransform range(vector<GTransform>& trans)
{
	GTransform result;
	double minX = MAXINT; double maxX = -MAXINT;
	double minY = MAXINT; double maxY = -MAXINT;

	for (int i = 0; i < trans.size(); i++)
	{
		GTransform& t = trans[i];
		minX = min(minX, t.dx);
		maxY = max(maxY, t.dy);
		maxY = max(maxY, t.dy);
	}

	result.dx = maxX - minX;
	result.dy = maxY - minY;

	return result;
}


double
muScore(vector<GTransform>& trans, double Dx, double Dy, double C2)
{
	// % Calculate differences between angles

	int sz = (int)trans.size();
	MatrixXd difAng(sz, sz);

	for (int c1 = 0; c1 < sz; c1++)
		for (int c2 = c1; c2 < sz; c2++)
		{
			difAng(c1, c2) = fmod(trans[c1].theta - trans[c2].theta + 2 * PI, 2 * PI);
			if (difAng(c1, c2) > PI)
				difAng(c1, c2) = 2 * PI - difAng(c1, c2);
		}

	GTransform temp = range(trans);

	double mu = 1 - C2 * (difAng.maxCoeff() / PI + temp.dx / Dx + temp.dy / Dy);

	return mu < 0 ? 0 : mu;
}


void CalcMeanTranslation(CFit& fit)
{
	int sz = (int)fit.m_ArcTrans.size();
	double tx = 0.0, ty = 0.0;
	for (int i = 0; i < sz; i++)
	{
		GTransform& trans = fit.m_ArcTrans[i];
		tx += trans.dx;
		ty += trans.dy;
	}
	fit.m_gFit.dx = tx / sz;
	fit.m_gFit.dy = ty / sz;

}


void
CalcAverageTheta(CFit& fit)
{
	double cang = 0.0, sang = 0.0;
	for (int c4 = 0; c4 < fit.m_Size; c4++)
	{
		cang += cos(fit.m_ArcTrans[c4].theta);
		sang += sin(fit.m_ArcTrans[c4].theta);
	}

	fit.m_gFit.theta = atan2(sang, cang);
	fit.m_gFit.dx = 0.0;
	fit.m_gFit.dy = 0.0;
}


void
CalcTranslation(CFit& fit, vector<CPlacement>& placements, vector<CTracker>& tracker)
{
	for (int c3 = 0; c3 < fit.m_Arcs.rows(); c3++)
	{
		const Curve& curve1 = Pieces[fit.m_Pieces(0)].m_Arcs(fit.m_Arcs(c3, 0)).m_Contour;
		Curve curve2 = TransformCurve(
			Pieces[fit.m_Pieces(1)].m_Arcs(fit.m_Arcs(c3, 1)).m_Contour,
			placements[tracker.back().m_Pc2Place[fit.m_Pieces(1)]].m_gLock);
		curve2.x.reverseInPlace();
		curve2.y.reverseInPlace();

		const Curve& sig1 = Pieces[fit.m_Pieces(0)].m_Arcs(fit.m_Arcs(c3, 0)).m_Signature;
		Curve sig2 = Orient_Reverse(Pieces[fit.m_Pieces(1)].m_Arcs(fit.m_Arcs(c3, 1)).m_Signature);

		Rigid_Motion_Translation(fit.m_ArcTrans[c3],
			curve1,
			curve2,
			sig1,
			sig2, fit.m_gFit.theta, params.beta);
	}
	fit.m_Score = muScore(fit.m_ArcTrans, Dx, Dy, params.C2);
	CalcMeanTranslation(fit);
}

void PlacePieces()
{
	size_t nPieces = Pieces.size();

	CPScore PScores(nPieces);
	vector<CPlacement> Placements;
	vector<CTracker> Tracker;
	Tracker.resize(1);

	//	Tracker.back().SetSize(nPieces);

	for (int i = 0; i < nPieces; i++)
	{
		Pieces[i].m_bActive = false;
	}

	int piece1;							// Start with the "heaviest" piece.
	double dummy = Weights.maxCoeff(&piece1);
//	piece1 = 11;									// ***BUGBUG*** We compute slightly different weights than Hoff.
	Pieces[piece1].m_bActive = true;

	Tracker.back().SetSize(nPieces);

	Tracker.back().m_PlacedPieces.resize(1);
	Tracker.back().m_PlacedPieces[0] = piece1;

	for (int i = 0; i < nPieces; i++)
	{
		Tracker.back().m_RemainingPieces[i] = i;

		Tracker.back().m_ActiveArcs[i].resize(Pieces[i].m_Arcs.size());
		for (int j = 0; j < Pieces[i].m_Arcs.size(); j++)
			Tracker.back().m_ActiveArcs[i][j] = j;

		Tracker.back().m_ActivePoints[i].resize(Pieces[i].m_Contour.size());
		for (int j = 0; j < Pieces[i].m_Contour.size(); j++)
			Tracker.back().m_ActivePoints[i][j] = j;

		Tracker.back().m_Pc2Place[i] = (i == piece1) - 1;
	}

	Tracker.back().m_RemainingPieces.erase(Tracker.back().m_RemainingPieces.begin() + piece1);

	Placements.emplace_back(piece1, Eigen::Vector4d(), GTransform(), CFit());
	int c1 = (int)Placements.size();
	int j = 0;																// Parameter sequence index

// c1	counts # of placed pieces
// c2	counts unplaced pieces.  Index into tracker.m_RemainingPieces
// c3	counts placed pieces.  Index into tracker.m_PlacedPieces

	while (c1 < nPieces)
	{
		vector<CFit> Fits;

		for (int c2 = 0; c2 < (nPieces - c1); c2++)
		{
			for (int c3 = 0; c3 < c1; c3++)
			{
				//% --------------------------------------------------------------------------
				//% Find P-Scores between arcs of remaining pieces and arcs of placed pieces
				//%--------------------------------------------------------------------------

				int nPiece1 = Tracker.back().m_RemainingPieces[c2];
				int nPiece2 = Tracker.back().m_PlacedPieces[c3];
				if (PScores.IsEmpty(nPiece1, nPiece2))
				{
					PScores(nPiece1, nPiece2) =
						MatrixXd::Zero(Pieces[nPiece1].m_Arcs.size(), Pieces[nPiece2].m_Arcs.size());

					// c4	counts over active arcs in Piece1
					// c5	counts over active arcs in Piece2

					for (int c4 = 0; c4 < Tracker.back().m_ActiveArcs[nPiece1].size(); c4++)
					{
						int nArc1 = Tracker.back().m_ActiveArcs[nPiece1][c4];

						for (int c5 = 0; c5 < Tracker.back().m_ActiveArcs[nPiece2].size(); c5++)
						{
							int nArc2 = Tracker.back().m_ActiveArcs[nPiece2][c5];

							PScores(nPiece1, nPiece2)(nArc1, nArc2) =
								SignatureSimilarity(Pieces[nPiece1].m_Arcs(nArc1).m_Signature,
									 Orient_Reverse(Pieces[nPiece2].m_Arcs(nArc2).m_Signature),
									Dkappa);
						}
					}
				}
//				PScores.Display(nPiece1, nPiece2);
//				PScores.Display();

				//%--------------------------------------------------------------------------
				//% Find sequences of consecutive high P-Scores

				//%--------------------------------------------------------------------------
				MatrixXi included = MatrixXi::Zero(PScores(nPiece1, nPiece2).rows(), PScores(nPiece1, nPiece2).cols());
				vector<Vector2i> tArcs;

				int c4Max = (int) Tracker.back().m_ActiveArcs[nPiece1].size();

				for (int c4 = 0; c4 < c4Max; c4++)
				{
					int nArc1 = Tracker.back().m_ActiveArcs[nPiece1][c4];

					int c5Max = (int) Tracker.back().m_ActiveArcs[nPiece2].size();

					for (int c5 = 0; c5 < c5Max; c5++)
					{
						int nArc2 = Tracker.back().m_ActiveArcs[nPiece2][c5];

						if (!included(nArc1, nArc2)
							&& PScores(nPiece1, nPiece2)(nArc1, nArc2) >= params.seq[j].p0)
						{
							included(nArc1, nArc2) = 1;

							tArcs.resize(0);
							tArcs.push_back(Vector2i(nArc1, nArc2));

							VectorXi& aarcs1 = Tracker.back().m_ActiveArcs[nPiece1];
							VectorXi& aarcs2 = Tracker.back().m_ActiveArcs[nPiece2];
							size_t sz1 = aarcs1.size();
							size_t sz2 = aarcs2.size();

							// Look c6 elements ahead in aarcs1 and c6 elements behind in aarcs2.

							int c6 = 1;													// Offset, not an index

							while (PScores(nPiece1, nPiece2)(aarcs1((c4 + c6) % sz1), aarcs2((c5 - c6 + sz2) % sz2)) >= params.seq[j].p0)
							{
								if (abs(aarcs1((c4 + c6) % sz1) - aarcs1((c4 + c6 - 1) % sz1)) > 1
									|| abs(aarcs2((c5 - c6 + sz2) % sz2) - aarcs2((c5 - c6 + 1 + sz2) % sz2)) > 1)
									break;
								else
								{
									int a1 = (aarcs1(c4) + c6) % included.rows();
									int a2 = (int)((aarcs2(c5) - c6 + included.cols()) % included.cols());

									tArcs.push_back(Vector2i(a1, a2));
									included(a1, a2) = 1;

									c6++;
								}
							}

							// Look c7 elements behind in aarcs1 and c6 elements ahead in aarcs2.

							int c7 = 1;

							while (PScores(nPiece1, nPiece2)(aarcs1((c4 - c7 + sz1) % sz1), aarcs2((c5 + c7) % sz2)) >= params.seq[j].p0)
							{
								if (abs(aarcs1((c4 - c7 + sz1) % sz1) - aarcs1((c4 - c7 + 1 + sz1) % sz1)) > 1
									|| abs(aarcs2((c5 + c7) % sz2) - aarcs2((c5 + c7 - 1) % sz2)) > 1)
									break;
								else
								{
									int a1 = (int)((aarcs1(c4) - c7 + included.rows()) % included.rows());
									int a2 = (int)((aarcs2(c5) + c7) % included.cols());

									tArcs.insert(tArcs.begin(), Vector2i(a1, a2));
									included(a1, a2) = 1;

									c7++;
								}
							}
							int minIncluded = (int) min(included.rows(), included.cols());
							if ((c6 + c7 - 1 - minIncluded) > 0)
							{
								for (int c8 = 0; c8 < c6 + c7 - 1 - minIncluded; c8++)
								{
									CFit fit;
									fit.m_Size = minIncluded;
									fit.m_Pieces = Vector2i(nPiece1, nPiece2);
									Map<MatrixXi> t((int*)tArcs.data(), tArcs.size(), 2);
									fit.m_Arcs = t.block(c8, 0, minIncluded, 2);
									fit.m_Slot = 1;
									fit.m_Score = -1.0;		// Not yet computed
									Fits.push_back(fit);
								}
							}
							else
							{
								CFit fit;
								fit.m_Size = c6 + c7 - 1;
								fit.m_Pieces = Vector2i(nPiece1, nPiece2);
								Map<MatrixXi> t((int*)tArcs.data(), tArcs.size(), 2);
								fit.m_Arcs = t;
								fit.m_Slot = 1;
								fit.m_Score = -1.0;		// Not yet computed
								Fits.push_back(fit);
							}
							if (Fits.back().m_Size < params.seq[j].m0)
								Fits.pop_back();
						}
					}
				}
			}
		}
		
		qsort(Fits.data(), Fits.size(), sizeof(CFit), [](const void* p1, const void* p2)->int
		{
			CFit* pFit1 = (CFit*)p1;
			CFit* pFit2 = (CFit*)p2;
			return (pFit1->m_Size - pFit2->m_Size);
		});

		int c2 = 0;

		while (c2 < Fits.size())
		{
			CFit& fitc2 = Fits[c2];
			Vector2i& piecesc2 = fitc2.m_Pieces;

			if (Pieces[piecesc2(1)].m_bActive &&
				IsMember(piecesc2(0), Tracker.back().m_RemainingPieces) &&
				!AnyMatch(fitc2.m_Arcs.col(1), Tracker.back().m_InactiveArcs[piecesc2(1)]))
			{
			//% --------------------------------------------------------------------------
			//% Calculate transformation and mu - score
			//% --------------------------------------------------------------------------
				if (fitc2.m_Score < 0)						// No score, yet
				{
					// ***BUGBUG***  or maybe not.  Transforms are always zero'ed when created

					//fitc2.m_ArcTrans = MatrixX3d::Zero(fitc2.m_Arcs.rows(), 3);

					Vector2d Theta(0.0, 0.0);

					fitc2.m_gFit.dx = 0; fitc2.m_gFit.dy = 0; fitc2.m_gFit.theta = 0;
					fitc2.m_ArcTrans.resize(fitc2.m_Arcs.rows());

					for (int c3 = 0; c3 < fitc2.m_Arcs.rows(); c3++)
					{
						const Curve& curve1 = Pieces[piecesc2(0)].m_Arcs(fitc2.m_Arcs(c3, 0)).m_Contour;

						Curve curve2 = TransformCurve(
							Pieces[piecesc2(1)].m_Arcs(fitc2.m_Arcs(c3, 1)).m_Contour,
							Placements[Tracker.back().m_Pc2Place[piecesc2(1)]].m_gLock);
						curve2.x.reverseInPlace();
						curve2.y.reverseInPlace();

						const Curve& sig1 = Pieces[piecesc2(0)].m_Arcs(fitc2.m_Arcs(c3, 0)).m_Signature;
						Curve sig2 = Orient_Reverse(Pieces[piecesc2(1)].m_Arcs(fitc2.m_Arcs(c3, 1)).m_Signature);

						fitc2.m_ArcTrans[c3].theta = Rigid_Motion_Angle(curve1, curve2, sig1, sig2, params.beta);

						Theta(0) += cos(fitc2.m_ArcTrans[c3].theta);
						Theta(1) += sin(fitc2.m_ArcTrans[c3].theta);
					}

					fitc2.m_gFit.theta = std::atan2(Theta(1), Theta(0));

					// %Calculate translations with after rotation through average angle theta

					CalcTranslation(fitc2, Placements, Tracker);
				}

				if (fitc2.m_Score >= 0 && fitc2.m_Score < params.seq[j].mu0)
				{
					//% In this case, it may be that a subset of this collection
					//% of arc pairings will have a high enough mu - score.
					//% We need only check this possibility if trimming an arc
					//% will not drop the number of arc pairs below m_0

					if (fitc2.m_Size > params.seq[j].m0)
					{
						//% Create a temporary structure

						vector<CFit> tFits;
						tFits.resize(2);

						//% We generate smaller sequences of matched pairings by
						//% deleting the first and last arcs of Fits(c2).Note,
						//% however, that we only consider the result of deleting
						//% the last arc if Fits(c2) did not arise from deleting
						//% the first arc of a larger fit(i.e., if Fits(c2).Slot
						//% == 1).This avoids generating duplicate fits.

						//% Delete last arc

						if (fitc2.m_Slot)
						{
							CFit& fit = tFits[0];

							fit.m_Pieces = fitc2.m_Pieces;

							fit.m_Size = fitc2.m_Size - 1;
							fit.m_Arcs = fitc2.m_Arcs.topRows(fit.m_Size);
							fit.m_ArcTrans = fitc2.m_ArcTrans;
							fit.m_ArcTrans.pop_back();

							CalcAverageTheta(fit);
							CalcTranslation(fit, Placements, Tracker);
						}

						//% Delete first arc

						{
							CFit& fit = tFits[1];

							fit.m_Pieces = fitc2.m_Pieces;
							fit.m_Size = fitc2.m_Size - 1;
							fit.m_Arcs = fitc2.m_Arcs.bottomRows(fit.m_Size);
							fit.m_ArcTrans = fitc2.m_ArcTrans;
							fit.m_ArcTrans.erase(fit.m_ArcTrans.begin());

							CalcAverageTheta(fit);
							CalcTranslation(fit, Placements, Tracker);
							fit.m_Slot = 0;
						}

						//% The new fits in tFits are shorter, so we must figure
						//% out where to place them within Fits in order to
						//% preserve the sorting on that structure

						int c3 = c2 + 1;
						while (c3 < Fits.size() && Fits[c3].m_Size == Fits[c2].m_Size)
							c3++;

						// % Insert tFits into Fits.

//						if (c3 < Fits.size())
						Fits.insert(Fits.begin() + c3, tFits.begin(), tFits.end());
					}

					//% Delete the fit that has now been subdivided.

					Fits.erase(Fits.begin()+c2);

					//% Update counter variable

					c2--;
				}
				else
				{
					Lock(fitc2.m_gFit, Pieces[fitc2.m_Pieces(0,0)].m_Contour, Tracker.back().m_SolvedPuzzleBoundary, j);
				}

/*



*/

			}
		}
	}
}


//for (int c3 = 0; c3 < fitc2.m_Arcs.rows(); c3++)
//{
//	const Curve& curve1 = Pieces[piecesc2(0)].m_Arcs(fitc2.m_Arcs(c3, 0)).m_Contour;
//	Curve curve2 = TransformCurve(
//		Pieces[piecesc2(1)].m_Arcs(fitc2.m_Arcs(c3, 1)).m_Contour,
//		Placements[Tracker.back().m_Pc2Place[piecesc2(1)]].m_gLock);
//	curve2.x.reverseInPlace();
//	curve2.y.reverseInPlace();

//	const Curve& sig1 = Pieces[piecesc2(0)].m_Arcs(fitc2.m_Arcs(c3, 0)).m_Signature;
//	Curve sig2 = Orient_Reverse(Pieces[piecesc2(1)].m_Arcs(fitc2.m_Arcs(c3, 1)).m_Signature);

//	Rigid_Motion_Translation(fitc2.m_ArcTrans[c3],
//		curve1,
//		curve2,
//		sig1,
//		sig2, fitc2.m_gFit.theta, params.beta);
//}

//fitc2.m_Score = muScore(fitc2.m_ArcTrans, Dx, Dy, params.C2);
//CalcMeanTranslation(fitc2);



						//for (int c4 = 0; c4 < fit.m_Size; c4++)
						//{
						//	const Curve& curve1 = Pieces[fit.m_Pieces(0)].m_Arcs(fit.m_Arcs(c4, 0)).m_Contour;
						//	Curve curve2 = TransformCurve(
						//		Pieces[fit.m_Pieces[1]].m_Arcs(fit.m_Arcs(c4, 1)).m_Contour,
						//		Placements[Tracker.back().m_Pc2Place[fit.m_Pieces[1]]].m_gLock);
						//	curve2.x.reverseInPlace();
						//	curve2.y.reverseInPlace();

						//	const Curve& sig1 = Pieces[fit.m_Pieces[0]].m_Arcs(fit.m_Arcs(c4, 0)).m_Signature;
						//	Curve sig2 = Orient_Reverse(Pieces[fit.m_Pieces[1]].m_Arcs(fit.m_Arcs(c4, 1)).m_Signature);

						//	Rigid_Motion_Translation(fit.m_ArcTrans[c4],
						//		curve1,
						//		curve2,
						//		sig1,
						//		sig2, fit.m_gFit.theta, params.beta);
						//}

						//fit.m_Score = muScore(fit.m_ArcTrans, Dx, Dy, params.C2);
						//CalcMeanTranslation(fit);
