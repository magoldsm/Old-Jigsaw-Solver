#pragma once

#include <iostream>
#include <atomic>

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/types_c.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "Eigen/Dense"


#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/tick_count.h"




using Curve = Eigen::Matrix<double, -1, 2>;

#define	PI			3.14159265358




#define USE_TBB

#ifdef USE_TBB
#define FOR_START(vbl, min, max)	tbb::parallel_for(min, max, [&](int vbl) {
#define	FOR_END						} );
#else
#define FOR_START(vbl, min, max)	for (int vbl = min; vbl < max; vbl++) {
#define FOR_END						}
#endif

void Solver();
