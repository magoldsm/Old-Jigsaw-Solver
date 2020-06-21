#include "pch.h"
#include "CPiece.h"
#include "SignatureSimilarity.h"
#include "CParameters.h"
#include "Utilities.h"
#include "Jigsaw Solver.h"
#include "Lock.h"
#include "CProgress.h"

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
	int sz1 = (int) sig1.rows();
	int sz2 = (int) sig2.rows();

	VectorXd h = VectorXd::Zero(sz1);

	// Equations 3.4, 3.5

	// Calculate strengths of correspondence
	
	for (int i = 0; i < sz1; i++)
	{
		Matrix<double, 1, 2> t = sig1.row(i);
		VectorXd d = (sig2.rowwise() - t).array().square().rowwise().sum().sqrt();
		//VectorXd d = ((sig2.col(0) - t.x).array().square() + (sig2.col(1)-t.y).array().square()).sqrt();

		// 	Calculate and sum h over j

		for (int j = 0; j < sz2; j++)
		{
			if (d[j] < D)
			{
				h[i] += 1.0 / (pow((d[j] / (D - d[j])), pParams->m_nGamma) + pParams->m_dEpsilon);
			}
		}
	}

	// Equation 3.8

	VectorXd pPoint = h.array() / (h.array() + pParams->m_nC1);

	// Equation 3.10

	VectorXd kPow = sig1.col(0).array().pow(pParams->m_nAlpha);
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
	result.resizeLike(points);

	double c = cos(trans.theta);
	double s = sin(trans.theta);

	result.col(0) = (points.col(0)*c - points.col(1)*s).array() + trans.dx;
	result.col(1) = (points.col(0)*s + points.col(1)*c).array() + trans.dy;

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
	dx1.resizeLike(curve1.col(0));
	dx2.resizeLike(curve2.col(0));
	dy1.resizeLike(curve1.col(1));
	dy2.resizeLike(curve2.col(1));

	int sz1 = (int)dx1.rows();
	int sz2 = (int)dx2.rows();

	dx1.block(1, 0, sz1 - 2, 1) = curve1.col(0).topRows(sz1 - 2) - curve1.col(0).bottomRows(sz1 - 2);
	dx2.block(1, 0, sz2 - 2, 1) = curve2.col(0).topRows(sz2 - 2) - curve2.col(0).bottomRows(sz2 - 2);
	dy1.block(1, 0, sz1 - 2, 1) = curve1.col(1).topRows(sz1 - 2) - curve1.col(1).bottomRows(sz1 - 2);
	dy2.block(1, 0, sz2 - 2, 1) = curve2.col(1).topRows(sz2 - 2) - curve2.col(1).bottomRows(sz2 - 2);

	dx1[0] = dx1[1]; dx1[sz1 - 1] = dx1[sz1 - 2];
	dx2[0] = dx2[1]; dx2[sz2 - 1] = dx2[sz2 - 2];
	dy1[0] = dy1[1]; dy1[sz1 - 1] = dy1[sz1 - 2];
	dy2[0] = dy2[1]; dy2[sz2 - 1] = dy2[sz2 - 2];

	dx1.array() *= signature1.col(1).array();
	dy1.array() *= signature1.col(1).array();

	dx2.array() *= signature2.col(1).array();
	dy2.array() *= signature2.col(1).array();

	VectorXd ks1pow = signature1.col(1).array().pow(beta);
	VectorXd ks2pow = signature2.col(1).array().pow(beta);
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
	VectorXd ks1pow = signature1.col(1).array().pow(beta);
	VectorXd ks2pow = signature2.col(1).array().pow(beta);
	double ks1powsum = ks1pow.sum();
	double ks2powsum = ks2pow.sum();

	// Calculate center of mass weighted by kappa_s

	double cm1x = (curve1.col(0).array() * ks1pow.array()).sum() / ks1powsum;
	double cm1y = (curve1.col(1).array() * ks1pow.array()).sum() / ks1powsum;
	double cm2x = (curve2.col(0).array() * ks2pow.array()).sum() / ks2powsum;
	double cm2y = (curve2.col(1).array() * ks2pow.array()).sum() / ks2powsum;

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
		curve2.colwise().reverseInPlace();

		const Curve& sig1 = Pieces[fit.m_Pieces(0)].m_Arcs(fit.m_Arcs(c3, 0)).m_Signature;
		Curve sig2 = Orient_Reverse(Pieces[fit.m_Pieces(1)].m_Arcs(fit.m_Arcs(c3, 1)).m_Signature);

		Rigid_Motion_Translation(fit.m_ArcTrans[c3],
			curve1,
			curve2,
			sig1,
			sig2, fit.m_gFit.theta, pParams->m_nBeta);
	}
	fit.m_Score = muScore(fit.m_ArcTrans, Dx, Dy, pParams->m_nC2);
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
		Pieces[i].m_nActive = 0;
	}

	int piece1;							// Start with the "heaviest" piece.
	double dummy = Weights.maxCoeff(&piece1);
//	piece1 = 11;									// ***BUGBUG*** We compute slightly different weights than Hoff.
	Pieces[piece1].m_nActive = 1;

	Tracker.back().SetSize(nPieces);

	Tracker.back().m_PlacedPieces.resize(1);
	Tracker.back().m_PlacedPieces[0] = piece1;

	Tracker.back().m_SolvedPuzzleBoundary = Pieces[piece1].m_Contour;	// Initially, the solved puzzle is just piece1's contour

	Tracker.back().m_SPB_Pt2PcPt.resize(Pieces[piece1].m_Contour.rows(), 2);
	Tracker.back().m_SPB_Pt2PcPt.col(0).setConstant(piece1);
	for (int i = 0; i < Pieces[piece1].m_Contour.rows(); i++)
		Tracker.back().m_SPB_Pt2PcPt(i, 1) = i;


	for (int i = 0; i < nPieces; i++)
	{
		Tracker.back().m_RemainingPieces[i] = i;

		Tracker.back().m_ActiveArcs[i].resize(Pieces[i].m_Arcs.size());
		for (int j = 0; j < Pieces[i].m_Arcs.size(); j++)
			Tracker.back().m_ActiveArcs[i][j] = j;

		Tracker.back().m_ActivePoints[i].resize(Pieces[i].m_Contour.rows());
		for (int j = 0; j < Pieces[i].m_Contour.rows(); j++)
			Tracker.back().m_ActivePoints[i][j] = j;

		Tracker.back().m_Pc2Place[i] = (i == piece1) - 1;
	}

	Tracker.back().m_RemainingPieces.erase(Tracker.back().m_RemainingPieces.begin() + piece1);

	Placements.emplace_back(piece1, Eigen::Vector4d(), GTransform(), CFit());
	int c1 = (int)Placements.size();
	int nRPc = (int) nPieces - c1 + 1;

	int j = 0;																// Parameter sequence index

// c1	counts # of placed pieces
// c2	counts unplaced pieces.  Index into tracker.m_RemainingPieces
// c3	counts placed pieces.  Index into tracker.m_PlacedPieces

	Progress.RestartReport(PROGRESS_PLACING, true);

	while (c1 < nPieces)
	{
		vector<CFit> Fits;

		Progress.RestartReport(PROGRESS_COMPARING, true);

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
							&& PScores(nPiece1, nPiece2)(nArc1, nArc2) >= pParams->m_P0[j])
						{
							included(nArc1, nArc2) = 1;

							tArcs.resize(0);
							tArcs.push_back(Vector2i(nArc1, nArc2));

							vector<int>& aarcs1 = Tracker.back().m_ActiveArcs[nPiece1];
							vector<int>& aarcs2 = Tracker.back().m_ActiveArcs[nPiece2];
							size_t sz1 = aarcs1.size();
							size_t sz2 = aarcs2.size();

							// Look c6 elements ahead in aarcs1 and c6 elements behind in aarcs2.

							int c6 = 1;													// Offset, not an index

							while (PScores(nPiece1, nPiece2)(aarcs1[(c4 + c6) % sz1], aarcs2[(c5 - c6 + sz2) % sz2]) >= pParams->m_P0[j])
							{
								if (abs(aarcs1[(c4 + c6) % sz1] - aarcs1[(c4 + c6 - 1) % sz1]) > 1
									|| abs(aarcs2[(c5 - c6 + sz2) % sz2] - aarcs2[(c5 - c6 + 1 + sz2) % sz2]) > 1)
									break;
								else
								{
									int a1 = (aarcs1[c4] + c6) % included.rows();
									int a2 = (int)((aarcs2[c5] - c6 + included.cols()) % included.cols());

									tArcs.push_back(Vector2i(a1, a2));
									included(a1, a2) = 1;

									c6++;
								}
							}

							// Look c7 elements behind in aarcs1 and c6 elements ahead in aarcs2.

							int c7 = 1;

							while (PScores(nPiece1, nPiece2)(aarcs1[(c4 - c7 + sz1) % sz1], aarcs2[(c5 + c7) % sz2]) >= pParams->m_P0[j])
							{
								if (abs(aarcs1[(c4 - c7 + sz1) % sz1] - aarcs1[(c4 - c7 + 1 + sz1) % sz1]) > 1
									|| abs(aarcs2[(c5 + c7) % sz2] - aarcs2[(c5 + c7 - 1) % sz2]) > 1)
									break;
								else
								{
									int a1 = (int)((aarcs1[c4] - c7 + included.rows()) % included.rows());
									int a2 = (int)((aarcs2[c5] + c7) % included.cols());

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
							if (Fits.back().m_Size < pParams->m_M0[j])
								Fits.pop_back();
						}
					}
				}
				Progress(PROGRESS_COMPARING).m_Percent = 1.0*(c1*c2 + c3) / (c1*(nPieces - c1));
				Progress.UpdateReport();
				
			}
		}
		
		qsort(Fits.data(), Fits.size(), sizeof(CFit), [](const void* p1, const void* p2)->int
		{
			CFit* pFit1 = (CFit*)p1;
			CFit* pFit2 = (CFit*)p2;
			return (pFit1->m_Size - pFit2->m_Size);
		});

		int c2 = 0;
		bool bProgress = false;

		Progress.RestartReport(PROGRESS_CHECKING, true);

		// Check fits

		while (c2 < Fits.size())
		{
			CFit& fitc2 = Fits[c2];
			Vector2i& piecesc2 = fitc2.m_Pieces;

//        if(pieces(Fits(c2).Pieces(1, 2)).Active && ismember(Fits(c2).Pieces(1, 1), tracker(end).RPc) &&
//			~any(ismember(Fits(c2).Arcs(:, 2), tracker(end).IArcs{Fits(c2).Pieces(1, 2)})))

			if (Pieces[piecesc2(1)].m_nActive &&
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
						curve2.colwise().reverseInPlace();

						const Curve& sig1 = Pieces[piecesc2(0)].m_Arcs(fitc2.m_Arcs(c3, 0)).m_Signature;
						Curve sig2 = Orient_Reverse(Pieces[piecesc2(1)].m_Arcs(fitc2.m_Arcs(c3, 1)).m_Signature);

						fitc2.m_ArcTrans[c3].theta = Rigid_Motion_Angle(curve1, curve2, sig1, sig2, pParams->m_nBeta);

						Theta(0) += cos(fitc2.m_ArcTrans[c3].theta);
						Theta(1) += sin(fitc2.m_ArcTrans[c3].theta);
					}

					fitc2.m_gFit.theta = std::atan2(Theta(1), Theta(0));

					// %Calculate translations with after rotation through average angle theta

					CalcTranslation(fitc2, Placements, Tracker);
				}

				if (fitc2.m_Score >= 0 && fitc2.m_Score < pParams->m_MU0[j])
				{
					//% In this case, it may be that a subset of this collection
					//% of arc pairings will have a high enough mu - score.
					//% We need only check this possibility if trimming an arc
					//% will not drop the number of arc pairs below m_0

					if (fitc2.m_Size > pParams->m_M0[j])
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

					Indices tPiecePtIcs_3, tPiecePtIcs;
					GTransform gLock;
					Lock(fitc2.m_gFit, Pieces[fitc2.m_Pieces(0,0)].m_Contour, Tracker.back().m_SolvedPuzzleBoundary, 
						gLock, pParams->m_K3[j], tPiecePtIcs_3, Tracker.back().m_SP_Bdry_PtIcs_3, tPiecePtIcs, Tracker.back().m_SP_Bdry_PtIcs);

					// Compute Scores

					double q1 = 1.0 * Tracker.back().m_SP_Bdry_PtIcs_3.size() / Tracker.back().m_SP_Bdry_PtIcs.size();
					double q2 = 0;
					int tn = (int) tPiecePtIcs_3.size();
					for (int c3 = 0; c3 < tn; c3++)
					{
						int c3p1 = (c3 + 1) % tn;
						int diff = abs(tPiecePtIcs_3[c3p1] - tPiecePtIcs_3[c3]);
						if (diff == 1 || diff == (Pieces[piecesc2[0]].m_Contour).size() - 1)
						{
							q2 += (Pieces[piecesc2[0]].m_Contour.row(tPiecePtIcs_3[c3p1]) - Pieces[piecesc2[0]].m_Contour.row(tPiecePtIcs_3[c3])).norm();
						}
					}
					q2 /= AverageLength;
					auto q3 = Gather(Pieces[piecesc2[0]].m_Signature, tPiecePtIcs_3).col(0).array().abs().sum()/(Dkappa*AverageSize);
						
					// Update progress report

					Progress(PROGRESS_CHECKING).m_Percent = 1.0*c2 / Fits.size();
					Progress.UpdateReport();

					// If scores are high enough, place the piece

					if ((pParams->m_dEta1*q1 + pParams->m_dEta2 * q3 > pParams->m_dQ1 || q2 > pParams->m_dQ2Star) && q2 > pParams->m_dQ2 && q3 > pParams->m_dQ3)
					{
						// Add new placement to placements variable
							
						Placements[c1].m_Score << q1, q2, q3, pParams->m_dEta1*q1 + pParams->m_dEta2 * q3;
						Placements[c1].m_nPiece = fitc2.m_Pieces[0];
						Placements[c1].m_gLock = gLock;
						Placements[c1].m_Fit = fitc2;

						// Update tracker variable

						Tracker.push_back(Tracker.back());
						Tracker.back().m_Pc2Place[Placements[c1].m_nPiece] = c1;

						// Mark arcs / points as used or remaining and introduce neighbors / activate

						vector<int>& inactive = Tracker.back().m_InactiveArcs[piecesc2[0]];

						for (int c3 = 0; c3 < tPiecePtIcs.size(); c3++)
						{
							auto arcno = Pieces[piecesc2[0]].m_Pt2Arc(tPiecePtIcs[c3]);
							if (arcno)
								inactive.push_back(arcno);
						}

						Placements[c1].m_Neighbors.resize(0);
						Pieces[fitc2.m_Pieces(0)].m_nActive = 2;

						vector<int>& neighbors = Placements[c1].m_Neighbors;

						for (int c3 = 0; c3 < Tracker.back().m_SP_Bdry_PtIcs.size(); c3++)
						{
							int SPIc3 = Tracker.back().m_SP_Bdry_PtIcs[c3];

							// Introduce new neighbors

							int newNeighbor = Tracker.back().m_SPB_Pt2PcPt(SPIc3, 0);

							if (!IsMember(newNeighbor, neighbors))
							{
								neighbors.push_back(newNeighbor);
								Pieces[newNeighbor].m_nActive = 2;

								vector<int>& toSetActive = Placements[Tracker[c1].m_Pc2Place[newNeighbor]].m_Neighbors;
								for (int i = 0; i < toSetActive.size(); i++)
									Pieces[toSetActive[i]].m_nActive = 2;

								toSetActive.push_back(fitc2.m_Pieces(0));

							}

							// Mark arcs as inactive

							auto arc = Pieces[newNeighbor].m_Pt2Arc(Tracker.back().m_SPB_Pt2PcPt(SPIc3, 1), 1);

							if (arc)
							{
								Tracker.back().m_InactiveArcs[newNeighbor].push_back(arc);
							}

							// Mark points as inactive

							Tracker.back().m_InactivePoints[newNeighbor].push_back(Tracker.back().m_SPB_Pt2PcPt(SPIc3, 1));

						}

						auto& fitc2piece = fitc2.m_Pieces(0);

						Tracker.back().m_InactivePoints[fitc2piece] = tPiecePtIcs;

						// Note: Lock returns tPiecePtIcs in sorted order.  Convenient!
						//
						// ActivePoints is the complement of the InactivePoints

						vector<int>& activePoints = Tracker.back().m_ActivePoints[fitc2piece];

						for (int i = 0; i < Pieces[fitc2piece].m_Contour.rows(); i++)
							activePoints[i] = i;
						for (int i = (int) tPiecePtIcs.size() - 1; i >= 0; i--)
							activePoints.erase(activePoints.begin() + tPiecePtIcs[i]);

						// Same with arcs, except might not be sorted and will have duplicates!

						vector<int>& activeArcs = Tracker.back().m_ActiveArcs[fitc2piece];

						AllExcept(activeArcs, Pieces[fitc2piece].m_Arcs.rows(), Tracker.back().m_InactiveArcs[fitc2piece]);

						//for (int i = 0; i < Pieces[fitc2piece].m_Arcs.rows(); i++)
						//	activeArcs[i] = i;
						//vector<int> inactiveArcs = Tracker.back().m_InactiveArcs[fitc2piece];
						//auto icmp = [](const void* p1, const void* p2) { return *(int*)p1 - *(int*)p2; };
						//qsort(inactiveArcs.data(), inactiveArcs.size(), sizeof(int), icmp);
						//auto ucmp = [](int i, int j) { return i == j; };
						//unique(inactiveArcs.begin(), inactiveArcs.end(), ucmp);
						//for (int i = inactiveArcs.size() - 1; i >= 0; i--)
						//	activeArcs.erase(activeArcs.begin() + inactiveArcs[i]);

						for (int c3 = 0; c3 < c1 - 1; c3++)
						{
							int placedPiece = Tracker.back().m_PlacedPieces[c3];
							vector<int>& activeArcs = Tracker.back().m_ActiveArcs[placedPiece];
							vector<int>& activePoints = Tracker.back().m_ActivePoints[placedPiece];

							// Calculate AArcs

							AllExcept(activeArcs, Pieces[placedPiece].m_Arcs.size(), Tracker.back().m_InactiveArcs[placedPiece]);

							// Calculate APts

							AllExcept(activePoints, Pieces[placedPiece].m_Contour.rows(), Tracker.back().m_InactivePoints[placedPiece]);
						}

						// Calculate new Fragment

						// Remove elements from the solved puzzle boundary

						Tracker.back().m_SolvedPuzzleBoundary = RemoveElements(Tracker.back().m_SolvedPuzzleBoundary, Tracker.back().m_SP_Bdry_PtIcs);

						Curve tPiece = Pieces[fitc2.m_Pieces(0)].m_Contour;
						vector<int> tTrack;
						AllExcept(tTrack, tPiece.size(), tPiecePtIcs);
						tPiece = RemoveElements(tPiece, tPiecePtIcs);

						Append(Tracker.back().m_SolvedPuzzleBoundary, tPiece);

						//Curve joined(Tracker.back().m_SolvedPuzzleBoundary.rows() + tPiece.rows(), 2);
						//joined << Tracker.back().m_SolvedPuzzleBoundary, TransformCurve(tPiece, gLock);
						//Tracker.back().m_SolvedPuzzleBoundary = joined;

						Tracker.back().m_SPB_Pt2PcPt = RemoveElements(Tracker.back().m_SPB_Pt2PcPt, Tracker.back().m_SP_Bdry_PtIcs);

						MatrixX2i tTrackPlus(tTrack.size(), 2);
						tTrackPlus.col(0).setConstant(Placements[c1].m_nPiece);

						for (int i = 0; i < tTrack.size(); i++)
						{
							tTrackPlus(i, 1) = tTrack[i];
						}

						Append(Tracker.back().m_SPB_Pt2PcPt, tTrackPlus);

						// Mark pieces as used or remaining

						Tracker.back().m_PlacedPieces.push_back(Placements[c1].m_nPiece);
						AllExcept(Tracker.back().m_RemainingPieces, nPieces, Tracker.back().m_PlacedPieces);

						Progress(PROGRESS_PLACING) = (int) (c1 - nPieces + nRPc) / nRPc;
						Progress.UpdateReport();

						bProgress = true;
						c1++;
					}
				}
			}
		}
		
		Progress.RestartReport(PROGRESS_CHECKING, false);
		Progress(PROGRESS_COMPARING) = 0;
		Progress.UpdateReport();

		// If the algorithm has dead - ended, increase depth in paraemter sequence
		// or terminate algorithm if this is not possible

		if (!bProgress)
		{
			if (j < pParams->m_nJStar)
			{
				j++;

				// Activate pieces

				for (int c4 = 0; c4 < Tracker.back().m_RemainingPieces.size(); c4++)
				{
					Pieces[Tracker.back().m_PlacedPieces[c4]].m_nActive = 1;
				}
			}
			else
			{
				// ***BUGBUG*** Dead End!  Also, plot!
				__debugbreak();
				return;
			}
		}
		else
		{
			j = 1;
			// Deactivate pieces
			for (int c4 = 0; c4 < Tracker.back().m_PlacedPieces.size(); c4++)
			{
				Pieces[Tracker.back().m_PlacedPieces[c4]].m_nActive = 
					max(0, Pieces[Tracker.back().m_PlacedPieces[c4]].m_nActive - 1);
			}
		}

		/*
		*/
	}
}


