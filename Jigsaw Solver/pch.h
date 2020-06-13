// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

#ifndef PCH_H
#define PCH_H

#include <iostream>

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/types_c.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "Eigen/Dense"


#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/tick_count.h"

#if 0
class Curve
{
public:
	Curve() : x(), y() {}
	Curve(int sz) : x(sz), y(sz) {}
	Curve(const Curve& c) : x(), y() {
		x = c.x;
		y = c.y;
	}
	~Curve() {}

	Curve& operator=(const Curve& c)
	{
		if (this == &c)
			return *this;
		x = c.x;
		y = c.y;
		return *this;
	}

	cv::Point2d operator[](int i) const
	{
		return cv::Point2d(x[i], y[i]);
	}

	size_t size() const { return x.size(); }
	void resize(size_t sz) 
	{ 
		x.resize(sz); 
		y.resize(sz);
	}
	union
	{
		Eigen::VectorXd		x;
		Eigen::VectorXd		kappa;
	};
	union
	{
		Eigen::VectorXd		y;
		Eigen::VectorXd		kappas;
	};
};
#else

using Curve = Eigen::Matrix<double, -1, 2>;

#endif


#define	PI			3.14159265358


#define USE_TBB

#ifdef USE_TBB
#define FOR_START(vbl, min, max)	tbb::parallel_for(min, max, [&](int vbl) {
#define	FOR_END						} );
#else
#define FOR_START(vbl, min, max)	for (int vbl = min; vbl < max; vbl++) {
#define FOR_END						}
#endif
#endif //PCH_H

/*
class Curve
{
public:
	Curve() : x(NULL, 0), y(NULL, 0) {}
	Curve(int sz) : x(NULL, 0), y(NULL, 0), curve(sz, sz)
	{
		Remap();
	}
	Curve(const Curve& c) : x(NULL, 0), y(NULL, 0) {
		curve = c.curve;
		Remap();
		//x = c.x;
		//y = c.y;
	}
	~Curve() {}

	Curve& operator=(const Curve& c)
	{
		if (this == &c)
			return *this;

		curve = c.curve;
		Remap();
		//x = c.x;
		//y = c.y;
		return *this;
	}

	cv::Point2d operator[](int i) const
	{
		return cv::Point2d(x[i], y[i]);
	}

	size_t size() const { return x.size(); }
	void resize(size_t sz)
	{
		curve.resize(sz, sz);
		Remap();
		//x.resize(sz);
		//y.resize(sz);
	}

	Eigen::Matrix2Xd					curve;

	union
	{
		Eigen::Map<Eigen::VectorXd>		x;
		Eigen::Map<Eigen::VectorXd>		kappa;
	};
	union
	{
		Eigen::Map<Eigen::VectorXd>		y;
		Eigen::Map<Eigen::VectorXd>		kappas;
	};

private:
	void Remap()
	{
		//new (&x) Eigen::Map<Eigen::VectorXd>(curve.row(0).data());
		//new (&y) Eigen::Map<Eigen::VectorXd>(curve.row(1).data());
	}
};
*/