// Jigsaw Solver.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>

#include "Savitsy-Golay.h"
#include "CParameters.h"
#include "Jigsaw Solver.h"
#include "Euclidean Signature.h"
#include "BivertexArcDecomposition.h"
#include "Utilities.h"
#include "SignatureSimilarity.h"
#include "CProgress.h"


using namespace std;
using namespace Eigen;
using namespace cv;

std::vector<CPiece> Pieces;
Eigen::VectorXd Weights;
VectorXd smoothVec, d1Vec, d2Vec;
double Dx, Dy, Dkappa, Dkappas;
double AverageLength;
long AverageSize;
CProgress Progress(4);


int main()
{
	LARGE_INTEGER liStart, liFrequency, liSG, liEuclid, liBVD, liTotal;

	QueryPerformanceFrequency(&liFrequency);

	//bool bSuccess = ReadPuzzle(Pieces, "Hoff - 1 piece.csv");
	bool bSuccess = ReadPuzzle(Pieces, "Hoff-2 pieces.csv");
	//CPiece piece = Pieces[35];
	//Pieces.resize(1);
	//Pieces[0] = piece;
	int nPieces = (int) Pieces.size();
	//double min, max;

	//MyNorm(Pieces[0].m_Contour.x, min, max);
	//MyNorm(Pieces[0].m_Contour.y, min, max);


	QueryPerformanceCounter(&liStart);

	GenerateSGVector(0, params.sgOrder, params.sgWindow, smoothVec);
	GenerateSGVector(1, params.sgOrder, params.sgWindow, d1Vec);
	GenerateSGVector(2, params.sgOrder, params.sgWindow, d2Vec);

	QueryPerformanceCounter(&liSG);
	cout << "Savitsky-Golay " << (liSG.QuadPart - liStart.QuadPart) / (1.0*liFrequency.QuadPart) << " seconds" << endl;

	/*
		Mat img = imread("Single Contour.jpg");

		Mat gray;
		cvtColor(img, gray, CV_BGR2GRAY);
		threshold(gray, gray, 128, 255, THRESH_BINARY_INV);
		//	imshow("grey", gray);

		vector<vector<Point> > Contours;
		findContours(gray, Contours, CV_RETR_EXTERNAL, CHAIN_APPROX_NONE);
		vector<Point>& contour = Contours[0];
		Rect rect = boundingRect(gray);
		int sz = (int)contour.size();

		// ***BUGBUG*** This code needs to be updated for multiple pieces.

		Pieces.resize(1);
		Pieces[0] = contour;
	*/

	AverageSize = 0;
	AverageLength = 0.0;
	Weights.resize(nPieces);
	atomic<int> nDone(0);
	Progress.RestartReport(PROGRESS_EUCLID, true);

	FOR_START(i, 0, nPieces)
		CalcEuclideanSignature(Pieces[i], smoothVec, d1Vec, d2Vec);
		nDone += 1;
		Progress[PROGRESS_EUCLID].m_Percent = 1.0*nDone / nPieces;
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

	AverageSize /= (long) Pieces.size();
	AverageLength /= Pieces.size();

	//for (int i = 0; i < Pieces[0].size(); i++)
	//{
	//	cout << Pieces[0].m_Contour.x[i] << "|" << Pieces[0].m_Contour.y[i] << "|" << Pieces[0].m_Signature.x[i] << "|" << Pieces[0].m_Signature.y[i] << endl;
	//}

	// Determine the characteristic size of the puzzle piece.
	// Basically, find a bounding rectangle that covers all piece contours and one that covers all their signatures

	MaxD(Pieces, &CPiece::m_Contour, Dx, Dy);
	MaxD(Pieces, &CPiece::m_Signature, Dkappa, Dkappas);

	// We use these to determine delta0 and delta1.  These are the amount of "noise" we tolerate
	// as we're dividing the pieces into their Bivertex Decompositions

	double delta0 = Dkappas / params.lambda0;
	double delta1 = Dkappa / params.lambda1;

	//vector<BATable> BA;
	//vector<Eigen::VectorXi> Pt2Arc;
	//BA.resize(Pieces.size());
	//Pt2Arc.resize(Pieces.size());

	nDone = 0;
	Progress.RestartReport(PROGRESS_BIVERTEX, true);

	FOR_START(i, 0, nPieces)
		BivertexArcDecomposition(Pieces[i], delta0, delta1/*, BA[i], Pt2Arc[i]*/);
		nDone++;
		Progress[PROGRESS_BIVERTEX].m_Percent = 1.0*nDone / nPieces;
		Progress.UpdateReport();
	FOR_END

	QueryPerformanceCounter(&liBVD);
	cout << "Bivertex Decomp " << (liBVD.QuadPart - liEuclid.QuadPart) / (1.0*liFrequency.QuadPart) << " seconds" << endl;

	Progress.RestartReport(PROGRESS_EUCLID, true);
	Progress[PROGRESS_EUCLID].m_Percent = 1.0;
	Progress.UpdateReport();
	Progress.RestartReport(PROGRESS_BIVERTEX, true);
	Progress[PROGRESS_BIVERTEX].m_Percent = 1.0;
	Progress.UpdateReport();

	//for (int i = 0; i < nPieces; i++)
	//{
	//	PlotDecomposition(Pieces[i], i);
	//	PlotKappasDecomposition(Pieces[i], delta0, i);
	//	waitKey(1);
	//}

	//waitKey(0);

	PlacePieces();
	QueryPerformanceCounter(&liTotal);
	cout << "Total Runtime " << (liTotal.QuadPart - liStart.QuadPart) / (1.0*liFrequency.QuadPart) << " seconds" << endl;
	waitKey(0);
}

