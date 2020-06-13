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


using namespace std;
using namespace Eigen;
using namespace cv;

std::vector<CPiece> Pieces;
Eigen::VectorXd Weights;
VectorXd smoothVec, d1Vec, d2Vec;
double Dx, Dy, Dkappa, Dkappas;

static VectorXd test(MatrixX2d x)
{
	VectorXd res(x.rows());
	
	res = x.col(0);

	return res;
}

int main()
{
	MatrixX2d x(5, 2);
	x << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

	cout << x << endl;

	VectorXd y = test(x);

	cout << y << endl;

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

	long AverageSize = 0;
	double AverageLength = 0.0;
	Weights.resize(nPieces);

	FOR_START(i, 0, nPieces)
		CalcEuclideanSignature(Pieces[i], smoothVec, d1Vec, d2Vec);
		//vector<Curve> test;
		//test.push_back(Pieces[i].m_Contour);
		//test.push_back(Pieces[i].m_Signature);

		//PlotContours(test, "test", false);
		//waitKey(0);
	FOR_END

	for (int i = 0; i < nPieces; i++)
	{
		Curve& c = Pieces[i].m_Contour;
		VectorXd xm1, ym1;
		circShift(c.x, xm1, -1);
		circShift(c.y, ym1, -1);

		Weights[i] = Pieces[i].m_Signature.kappa.array().abs().sum();

		AverageSize += (long)Pieces[i].m_Contour.size();
		AverageLength += ((c.x - xm1).array().square() + (c.y - ym1).array().square()).sqrt().sum();
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

	FOR_START(i, 0, nPieces)
		BivertexArcDecomposition(Pieces[i], delta0, delta1/*, BA[i], Pt2Arc[i]*/);
	FOR_END

	QueryPerformanceCounter(&liBVD);
	cout << "Bivertex Decomp " << (liBVD.QuadPart - liEuclid.QuadPart) / (1.0*liFrequency.QuadPart) << " seconds" << endl;

	for (int i = 0; i < nPieces; i++)
	{
		//PlotDecomposition(Pieces[i], i);
		//PlotKappasDecomposition(Pieces[i], delta0, i);
		waitKey(1);
	}

	waitKey(0);

	PlacePieces();
	QueryPerformanceCounter(&liTotal);
	cout << "Total Runtime " << (liTotal.QuadPart - liStart.QuadPart) / (1.0*liFrequency.QuadPart) << " seconds" << endl;
	waitKey(0);
}

