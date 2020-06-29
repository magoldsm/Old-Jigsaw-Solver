#include "pch.h"
#include "CProgress.h"
#include "CParameters.h"
#include "Savitsy-Golay.h"
#include "Euclidean Signature.h"
#include "BivertexArcDecomposition.h"
#include "Utilities.h"
#include "SignatureSimilarity.h"
#include "CPiece.h"

using namespace std;
using namespace Eigen;
using namespace cv;

extern std::vector<CPiece> Pieces;
extern Eigen::VectorXd Weights;
extern Eigen::VectorXd smoothVec, d1Vec, d2Vec;
extern double AverageLength;
extern long AverageSize;
extern double Dx, Dy, Dkappa, Dkappas;


void Solver()
{
	LARGE_INTEGER liStart, liFrequency, liSG, liEuclid, liBVD, liTotal;

	QueryPerformanceFrequency(&liFrequency);

	bool bSuccess = ReadPuzzle(Pieces, pParams->m_szPath);
	int nPieces = (int)Pieces.size();

	QueryPerformanceCounter(&liStart);

	GenerateSGVector(0, pParams->m_nSGOrder, pParams->m_nSGWindow, smoothVec);
	GenerateSGVector(1, pParams->m_nSGOrder, pParams->m_nSGWindow, d1Vec);
	GenerateSGVector(2, pParams->m_nSGOrder, pParams->m_nSGWindow, d2Vec);

	QueryPerformanceCounter(&liSG);
	cout << "Savitsky-Golay " << (liSG.QuadPart - liStart.QuadPart) / (1.0*liFrequency.QuadPart) << " seconds" << endl;

	AverageSize = 0;
	AverageLength = 0.0;
	Weights.resize(nPieces);
	atomic<int> nDone(0);
	Progress.RestartReport(PROGRESS_EUCLID, true);

	FOR_START(i, 0, nPieces)
		CalcEuclideanSignature(Pieces[i], smoothVec, d1Vec, d2Vec);
		nDone += 1;
		Progress(PROGRESS_EUCLID).m_Percent = 1.0*nDone / nPieces;
		Progress.UpdateReport();
	FOR_END

		for (int i = 0; i < nPieces; i++)
		{
			Curve& c = Pieces[i].m_Contour;
			Curve cm1;
			cm1.resizeLike(c);

			circShift(c, cm1, -1);

			Weights[i] = Pieces[i].m_Signature.col(0).array().abs().sum();

			AverageSize += (long)Pieces[i].m_Contour.size();
			AverageLength += ((c.col(0) - cm1.col(0)).array().square() + (c.col(1) - cm1.col(1)).array().square()).sqrt().sum();
		}

	QueryPerformanceCounter(&liEuclid);
	cout << "Eucliean Sigs " << (liEuclid.QuadPart - liSG.QuadPart) / (1.0*liFrequency.QuadPart) << " seconds" << endl;

	AverageSize /= (long)Pieces.size();
	AverageLength /= Pieces.size();


	// Determine the characteristic size of the puzzle piece.
	// Basically, find a bounding rectangle that covers all piece contours and one that covers all their signatures

	MaxD(Pieces, &CPiece::m_Contour, Dx, Dy);
	MaxD(Pieces, &CPiece::m_Signature, Dkappa, Dkappas);

	// We use these to determine delta0 and delta1.  These are the amount of "noise" we tolerate
	// as we're dividing the pieces into their Bivertex Decompositions

	double delta0 = Dkappas / pParams->m_nLambda0;
	double delta1 = Dkappa / pParams->m_nLambda1;

	nDone = 0;
	Progress.RestartReport(PROGRESS_BIVERTEX, true);

	FOR_START(i, 0, nPieces)
		BivertexArcDecomposition(Pieces[i], delta0, delta1/*, BA[i], Pt2Arc[i]*/);
	nDone++;
	Progress(PROGRESS_BIVERTEX).m_Percent = 1.0*nDone / nPieces;
	Progress.UpdateReport();
	FOR_END

	QueryPerformanceCounter(&liBVD);
	cout << "Bivertex Decomp " << (liBVD.QuadPart - liEuclid.QuadPart) / (1.0*liFrequency.QuadPart) << " seconds" << endl;

	Progress.RestartReport(PROGRESS_EUCLID, true);
	Progress(PROGRESS_EUCLID).m_Percent = 1.0;
	Progress.UpdateReport();
	Progress.RestartReport(PROGRESS_BIVERTEX, true);
	Progress(PROGRESS_BIVERTEX).m_Percent = 1.0;
	Progress.UpdateReport();

	if (pParams->m_bPlotBVD)
	{
		for (int i = 0; i < nPieces; i++)
		{
			PlotDecomposition(Pieces[i], i);
			PlotKappasDecomposition(Pieces[i], delta0, i);
			//	waitKey(1);
		}
	}

	waitKey(0);

	PlacePieces();
	QueryPerformanceCounter(&liTotal);
	cout << "Total Runtime " << (liTotal.QuadPart - liStart.QuadPart) / (1.0*liFrequency.QuadPart) << " seconds" << endl;
//	waitKey(0);

}