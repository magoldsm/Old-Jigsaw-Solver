#pragma once

using EuclideanSignature = Curve;

EuclideanSignature Orient_Reverse(EuclideanSignature& sig);

void MyNormm(Eigen::MatrixXd& mat, Eigen::Vector2d& dmin, Eigen::Vector2d& dmax);
void MyNormv(Eigen::VectorXd& vec, double& dmin, double& dmax);
void MyNormv(Eigen::VectorXd& vec);
void PlotContours(std::vector<Curve>& Curves, const char* windowName, bool bOverlay);
void Plot(const char* window, const Eigen::VectorXd& vec, double dDelta = 0.0, bool bAnimate = false);
void circShift(const Eigen::VectorXd& vin, Eigen::VectorXd& vout, int shift);
void circShift(const Curve& min, Curve& mout, int shift);
void MaxD(const std::vector<CPiece>& piece, Curve CPiece::*curve, double& width, double& height);
bool IsMember(int s, std::vector<int> v);
bool AnyMatch(const Eigen::VectorXi& v1, const std::vector<int>& v2);

Eigen::MatrixX2d Gather(const Eigen::MatrixX2d& src, std::vector<int> indices);
void AllExcept(std::vector<int>& dest, size_t size, std::vector<int> except);


template<class T>
T RemoveElements(const  T& src, std::vector<int> remove)
{
	auto icmp = [](const void* p1, const void* p2) { return *(int*)p1 - *(int*)p2; };
	auto ucmp = [](int i, int j) { return i == j; };

	qsort(remove.data(), remove.size(), sizeof(int), icmp);
	unique(remove.begin(), remove.end(), ucmp);

	T res;
	res.resize(src.rows() - remove.size(), 2);
	int ires = 0;
	int iremove = 0;

	for (int i = 0; i < src.rows(); i++)
	{
		if (i < remove[iremove])
			res(ires++) = src(i);

		else if (i == remove[iremove])
			iremove++;
	}
	return res;
}


template<class T>
inline void Copy(Eigen::Matrix<T, -1, 1>& dest, const std::vector<T>& src)
{
	dest.resize(src.size());
	memcpy(dest.data(), src.data(), sizeof(T)*src.size());
}

template<class T>
void Append(T& dst, const T&src)
{
	T joined(dst.rows() + src.rows(), dst.cols());
	joined << dst, src;
	dst = joined;
}
