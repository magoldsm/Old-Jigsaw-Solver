#pragma once

using EuclideanSignature = Curve;

EuclideanSignature Orient_Reverse(EuclideanSignature& sig);

void MyNorm(Eigen::MatrixXd mat, Eigen::Vector2d& dmin, Eigen::Vector2d& dmax);
void MyNorm(Eigen::VectorXd& vec, double& dmin, double& dmax);
void MyNorm(Eigen::VectorXd& vec);
void PlotContours(std::vector<Curve>& Curves, const char* windowName, bool bOverlay);
void Plot(const char* window, const Eigen::VectorXd& vec, double dDelta = 0.0, bool bAnimate = false);
void circShift(const Eigen::VectorXd& vin, Eigen::VectorXd& vout, int shift);
void MaxD(const std::vector<CPiece>& piece, Curve CPiece::*curve, double& width, double& height);
bool IsMember(int s, std::vector<int> v);
bool AnyMatch(const Eigen::VectorXi& v1, const Eigen::VectorXi& v2);

