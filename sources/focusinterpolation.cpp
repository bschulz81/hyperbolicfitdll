// Copyright(c) < 2021 > 
// <Benjamin Schulz> 
// Responsible for:
// The implementation of the library in C++, 
// the development of a parallelized RANSAC inspired algorithm that can removes outlier data in the function "focusposition_Regression",
// the implementation of a repeated median regression in the function "focusposition_Regression" instead of a simple regression,
// the implementation of various known outlier detection methods within the RANSAC 
// (MAD, S, Q, and T, biweight-midvariance estimators, Grubb's test for outliers, Peirce's criterion, percentage based midvariance)
// the implementation of the image ananysis algorithm developed by the user C.Y. Tan with the help of opencv and libfitsio libraries. The algorithm generates a power spectrum in fourier 
// space which can be used in the curve fitting library for determining the focus point.

// <Stephen King> 
// Responsible for:
// The Algorithm idea for the hyperbolic fit in the internal function "Regression" with the linear regression and least squares comparison, first suggested in
// https://aptforum.com/phpbb/viewtopic.php?p=25998#p25998 

// <Jim Hunt>
// Responsible for:
// Some testing of the library,
// the algorithm idea for "findbackslash_Regression", first suggested in
// https://aptforum.com/phpbb/viewtopic.php?p=26265#p26265 

// <C.Y. Tan>  
// Responsible for:
// The development of the algorithm that analyzes the image data via Fourier transform in the function image::fouriertransform. 
// Its output can be fed into a slightly modified curve fitting algorithm.
// The suggestion in https://aptforum.com/phpbb/viewtopic.php?p=26471#p26471 to throw out outlier data by comparison of the error with the Standard-deviation.
// (Note that this old idea is now supplemented by more advanced methods, since the average and standard deviation are not robust.)

// The library makes use of an imaging analysis algorithm that was developed by C. Y. Tan (main author) with some contributions from B. Schulz, which is yet to be published.

// The library makes use of an algorithm for student's distribution, which can be found at
// Smiley W. Cheng, James C. Fu, Statistics & Probability Letters 1 (1983), 223-227

// The library also makes use of Peirce's outlier test. An algorithm for this method was developed in
// Gould, B. A, Astronomical Journal, vol. 4, iss. 83, p. 81 - 87 (1855).

// The library also has the possibility to use MAD, S, Q and T estimators. 
// These estimators are extensively described in 

// Peter J. Rousseeuw, Christophe Croux, Alternatives to the Median-Absolute Deviation
// J. of the Amer. Statistical Assoc. (Theory and Methods), 88 (1993),p. 1273,

// Christophe Croux and Peter J.Rousseeuw, Time-effcient algorithms for two highly robust estimators of scale, 
// In: Dodge Y., Whittaker J. (eds) Computational Statistics. Physica, Heidelberg, https ://doi.org/10.1007/978-3-662-26811-7_58

// According to Rousseeuw and Croux, the S estimator has an advantange over the MAD estimator because it can also work for asymmetric distributions.
// The same authors note that the Q estimator is better optimized than the S estimator for small sample sizes.
// 
// The library also can make use of the biweight midvariance estimator that was described in 
// T. C. Beers,K. Flynn and K. Gebhardt,  Astron. J. 100 (1),32 (1990)

// Within the Ransac, the offers the possibility to fit the data using Siegel's repeated median slope from 
// Siegel, Andrew, Technical Report No. 172, Series 2 By Department of Statistics Princeton University: Robust Regression Using Repeated Medians (September 1980)

// In practice, repeated median regression is a rather slow fitting method. 
// If median regression is not used, the RANSAC will use a faster linear regression algorithm for the 
// hyperbolic fit which was provided by Stephen King at https://aptforum.com/phpbb/viewtopic.php?p=25998#p25998).

// See 
// https://aptforum.com/phpbb/viewtopic.php?p=26587#p26587 for the early development history of this software. 
//
// A recent version of this software can be downloaded at
//
// https://github.com/bschulz81/hyperbolicfitdll/ 

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this softwareand associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and /or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <vector>
#include <chrono>
#include <algorithm>
#include "fitsio.h"

#include "focusinterpolation.h"
#include <random>  
#include <valarray>
#include <atomic>
#include <unordered_set>
#include <complex>
#include <cmath>

using namespace cv;
#if __cplusplus == 201703L
#include <mutex>
#include <execution>
#else
#include <omp.h>
#endif

struct maybe_inliner
{
	size_t point;
	double error;
};


inline void stdeviation(valarray<double>* errs, double* stdev, double* average, size_t s);
inline bool regression(double* ferr, valarray<long>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression, valarray<bool>* usepoint);
inline bool findmodel(valarray<long>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, valarray<bool>* usedpoint, size_t usedpoints, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method, size_t pointnumber, size_t pointnumberhalf);
inline bool Search_min_error(valarray<long>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression, vector<double>* errs, valarray<bool>* indicestouse);
inline double median(valarray<double>* arr, size_t n, size_t nhalf);
inline double lowmedian(valarray<double>* arr, size_t n);
inline double factorial(size_t n);
inline double peirce(size_t pointnumber, size_t numberofoutliers, size_t fittingparameters);
inline size_t binomial(size_t n, size_t k);
inline double crit(double alpha, size_t N);
inline double t(double alpha, size_t nu);
inline double cot(double x);
inline double G(double y, size_t nu);
inline double H(double y, size_t nu);
inline double Q_estimator(valarray<double>* err, size_t s);
inline double S_estimator(valarray<double>* err, size_t s);
inline double MAD_estimator(valarray<double>* err, double* m, size_t s, size_t shalf);
inline double onestepbiweightmidvariance(valarray<double>* err, double* median, size_t s, size_t shalf);
inline double T_estimator(valarray<double>* err, size_t s);


//Computes the standard deviation stdev and the average *average from an array *errs, espects the size of the array s
inline void stdeviation(valarray<double>* errs, double* stdev, double* average, size_t s)
{
	double sum = 0, devi = 0, t = 0;
	for (size_t i = 0; i < s; i++)
	{
		sum += (*errs)[i];
	}
	*average = sum / (double)s;

	if (stdev != NULL)
	{
		for (size_t i = 0; i < s; i++)
		{
			t = ((*errs)[i] - *average);
			t *= t;
			devi += t;
		}
		*stdev = sqrt(devi / s);
	}
}

// Computes the median of an array. Expects the size of the array and half of the size.
inline double median(valarray<double>* arr, size_t n, size_t nhalf)
{

#if __cplusplus == 201703L
	nth_element(std::execution::par, std::begin(*arr), std::begin(*arr) + nhalf, std::begin(*arr) + n);
#else
	nth_element(std::begin(*arr), std::begin(*arr) + nhalf, std::begin(*arr) + n);
#endif


	double  med = (*arr)[nhalf];
	if (n % 2 == 0)
	{
#if __cplusplus == 201703L
		auto max_it = max_element(std::execution::par, std::begin(*arr), std::begin(*arr) + nhalf);
#else
		auto max_it = max_element(std::begin(*arr), std::begin(*arr) + nhalf);
#endif

		med = (*max_it + med) / 2.0;
	}
	return med;
}

// computes the lower median of an array. expects the size of the array.
inline double lowmedian(valarray<double>* arr, size_t n)
{
	size_t m = (size_t)(floor(((double)n + 1.0) / 2.0) - 1.0);

#if __cplusplus == 201703L
	nth_element(std::execution::par, std::begin(*arr), std::begin(*arr) + m, std::begin(*arr) + n);
#else
	nth_element(std::begin(*arr), std::begin(*arr) + m, std::begin(*arr) + n);
#endif
	return (double)(*arr)[m];
}


// computes the fit of an array of points (x, line_y) to a hyperbola.
//  The minimum should lie between minfocus and maxfocus and is found to be at *focpos. it also returns the intercept and slope parameters as defined in the header.
// it accepts a use_median_regression parameter which, if true, will lead to a robust fit with Siegel's median slope. It also wants a boolean array usepoints that can mask certain points to be excluded from the fit.
// if the use_median_regression parameter is false, a linear regression published by Stephen King at https://aptforum.com/phpbb/viewtopic.php?p=25998#p25998  is used.
inline bool regression(double* ferr, valarray<long>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression, valarray<bool>* usepoint)
{
	*ferr = DBL_MAX;

	valarray<long> x2 = (*x)[*usepoint];
	valarray<double> line_y2 = (*line_y)[*usepoint];
	size_t usedpoints = x2.size();

	valarray<double>line_x2(usedpoints);


	if (use_median_regression)
	{
		//use Siegel's median slope with a naive algorithm

		valarray<double> stacks2(usedpoints);
		valarray<double> stacks1(usedpoints - 1);

		size_t halfsize = usedpoints / 2;

		valarray<double> stacki1(usedpoints);

		//for each possible focus point run siegel's median slope to fit a hyperbola. return the hyperbola that has a minimum such that it has the smallest error.
		for (long h = minfocus; h <= maxfocus; h++)
		{
			for (size_t i = 0; i < usedpoints; i++)
			{
				double k = (double)(x2)[i] - (double)h;
				line_x2[i] = k * k;
			}
			for (size_t i = 0; i < usedpoints; i++)
			{
				size_t q = 0;
				for (size_t j = 0; j < usedpoints; j++)
				{
					if (i != j)
					{
						double t = line_x2[j] - line_x2[i];
						if (t != 0)
						{
							stacks1[q] = (line_y2[j] - line_y2[i]) / t;
							q++;
						}
					}
				}
				stacks2[i] = median(&stacks1, q, q / 2);
			}

			double thisslope = median(&stacks2, usedpoints, halfsize);

			for (size_t n = 0; n < usedpoints; n++)
			{
				stacki1[n] = line_y2[n] - thisslope * line_x2[n];
			}
			double thisintercept = median(&stacki1, usedpoints, halfsize);

			double thiserr = 0;
			for (size_t n = 0; n < usedpoints; n++)
			{
				double	k = fabs(thisslope * line_x2[n] + thisintercept) - line_y2[n];
				thiserr += k * k;
			}

			if (thiserr < *ferr)
			{
				*fintercept = thisintercept;
				*fslope = thisslope;
				*ferr = thiserr;
				*focpos = h;
			}
		}
	}
	else
	{
		//make a linear regression fit for the hyperbola

		double sumy = 0;
		for (size_t n = 0; n < usedpoints; n++)
		{
			sumy += line_y2[n];
		}

		//for each possible focus position, fit a hyperbola. return the hyperbola with the focus position that has the smallest error to the points 
		for (long h = minfocus; h <= maxfocus; h++)
		{
			double sumx = 0, sumxy = 0, sumxx = 0;
			for (size_t n = 0; n < usedpoints; n++)
			{
				double k = (double)(x2)[n] - (double)h;
				k = k * k;
				line_x2[n] = k;
				sumx += k;
				sumxy += k * line_y2[n];
				sumxx += k * k;
			}

			double xaverage = sumx / (double)usedpoints;
			double yaverage = sumy / (double)usedpoints;

			double t = sumxx - sumx * xaverage;
			if (t == 0)
			{
				return false;
			}

			double thisslope = (sumxy - sumx * yaverage) / t;
			double thisintercept = yaverage - thisslope * xaverage;

			double thiserr = 0;

			for (size_t n = 0; n < usedpoints; n++)
			{
				double	k = fabs(thisslope * line_x2[n] + thisintercept) - line_y2[n];
				thiserr += k * k;
			}
			if (thiserr < *ferr)
			{
				*fintercept = thisintercept;
				*fslope = thisslope;
				*ferr = thiserr;
				*focpos = h;
			}
		}
	}

	*ferr = *ferr / (double)usedpoints;
	return true;
}

//starts the regression function and ensures that if the minimum point of the hyperbola is close to minfocus or maxfocus, the regression is started again on an enlarged area at the side close to minfocus or maxfocus. The size of this area
//is determined by the scale parameter.
inline bool Search_min_error(valarray<long>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression, valarray<bool>* indicestouse)
{

	if (!regression(err, x, line_y, minfocus, maxfocus, focpos, fintercept, fslope, use_median_regression, indicestouse))
	{
		return false;
	}
	long dist = abs(*focpos - maxfocus);
	if ((dist <= 10) && (scale > 1))
	{

		const long middle = lround(((double)maxfocus + (double)minfocus) / 2.0);
		minfocus = maxfocus;
		maxfocus = (long)((double)middle + ((double)maxfocus - (double)middle) * fabs(scale));
		long thisfocpos;
		double thiserr, thisslope, thisintercept;
		if (!regression(&thiserr, x, line_y, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope, use_median_regression, indicestouse))
		{
			return false;
		}

		if ((thiserr < *err))
		{
			*err = thiserr;
			*fintercept = thisintercept;
			*fslope = thisslope;
			*focpos = thisfocpos;
		}
	}
	dist = abs(*focpos - minfocus);
	if ((dist <= 10) && (scale > 1))
	{
		const long middle = lround(((double)maxfocus + (double)minfocus) / 2.0);
		maxfocus = minfocus;
		minfocus = (long)((double)middle - ((double)middle - (double)minfocus) * fabs(scale));
		long thisfocpos;
		double thiserr, thisslope, thisintercept;
		if (!regression(&thiserr, x, line_y, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope, use_median_regression, indicestouse))
		{
			return false;
		}
		if ((thiserr < *err))
		{
			*err = thiserr;
			*fintercept = thisintercept;
			*fslope = thisslope;
			*focpos = thisfocpos;
		}
	}
	return true;
}

//computes the factorial
_inline double factorial(size_t n)
{
	double ret = 1.00;
	for (size_t i = 2; i <= n; ++i)
		ret *= (double)i;
	return ret;
}

//helper function for student's t distribution 
//from the algorithm in Smiley W. Cheng, James C. Fu, Statistics & Probability Letters 1 (1983), 223-227
_inline double H(double y, size_t nu)
{
	double sum = 0;

	for (size_t j = 1; j <= size_t((nu + 1) / 2 - 1); j++)
	{
		sum += factorial((size_t)j) * factorial(size_t(j - 1.0)) / (pow(4, -((double)j)) * factorial((size_t)2.0 * j)) * pow((1 + y * y / nu), -((double)j));
	}
	return y / (2 * sqrt((double)nu)) * sum;
}

//helper function for student's t distribution from the algorithm in
// Smiley W. Cheng, James C. Fu, Statistics & Probability Letters 1 (1983), 223-227
_inline double G(double y, size_t nu)
{
	double sum = 0;
	for (size_t j = 0; j <= nu / 2 - 1; j++)
	{
		sum += factorial((size_t)2.0 * j) / (pow(2.0, 2.0 * j) * pow(factorial(j), 2)) * pow(1 + y * y / nu, -((double)j));
	}
	return sum;
}
//computes the cotangent
inline double cot(double x)
{
	return  cos(x) / sin(x);
}


//computes the student t distribution from the algorithm in
// // Smiley W. Cheng, James C. Fu, Statistics & Probability Letters 1 (1983), 223-227
inline double t(double alpha, size_t nu)
{
	const double PI = 3.14159265358979323846;
	if (nu % 2 != 0)
	{
		double zeta;
		double zeta1 = cot(alpha * PI);
		do
		{
			zeta = zeta1;
			zeta1 = sqrt((double)nu) * cot(alpha * PI + H(zeta, nu));
		} while (fabs(zeta1 - zeta) > 0.0001);
		return zeta1;
	}
	else
	{
		double zeta;
		double zeta1 = sqrt(2.0 * pow(1.0 - 2.0 * alpha, 2.0) / (1 - pow(1.0 - 2.0 * alpha, 2.0)));
		do
		{
			zeta = zeta1;
			zeta1 = pow(1.0 / nu * (pow(G(zeta, nu) / (1.0 - 2.0 * alpha), 2.0) - 1.0), -0.5);
		} while (fabs(zeta1 - zeta) > 0.0001);
		return zeta1;
	}
}
//computes the critical values of the student t distribution
inline double crit(double alpha, size_t N)
{
	double temp = pow(t(alpha / (2.0 * (double)N), N - 2), 2.0);
	return ((double)N - 1.0) / sqrt((double)N) * sqrt(temp / ((double)N - 2.0 + temp));
}

//computes the peirce criterion for a given pointnumber with a given number of outliers and a number of fitting parameters
inline double peirce(size_t pointnumber, size_t numberofoutliers, size_t fittingparameters)
{
	double n = (double)numberofoutliers, N = (double)pointnumber, p = (double)fittingparameters;
	double xa = 0.0, R = 1.0, diff1 = N - n;
	double xb;
	if (N > 1)
	{
		double lnQ = n * log(n) / N - log(N) + diff1 * log(diff1) / N;
		do
		{
			double lambda = exp((N * lnQ - n * log(R)) / (diff1));
			xb = xa;
			xa = 1.0 + (diff1 - p) / n * (1.0 - lambda * lambda);
			if (xa < 0)
			{
				xa = 0.0;
			}
			else
			{
				R = exp((xa - 1) / 2.0) * erfc(sqrt(xa) / sqrt(2.0));
			}
		} while (fabs(xa - xb) > (1.0e-5));
	}
	return xa;
}

//computes the binomial coeffocient n choose k
inline size_t binomial(size_t n, size_t k)
{
	if (k == 0)
	{
		return 1;
	}
	else if (n == k)
	{
		return 1;
	}
	else
	{
		if (n - k < k)
		{
			k = n - k;
		}
		double prod1 = 1.0, prod2 = (double)n;
		for (size_t i = 2; i <= k; i++)
		{
			prod1 *= i;
			prod2 *= ((double)n + 1.0 - (double)i);
		}
		return(size_t)(prod2 / prod1);
	}
}

//computes the Q estimator of an array of values with a naive algorithm. The array should have a size of s
// The q estimator was published in Peter J. Rousseeuw, Christophe Croux, Alternatives to the Median-Absolute Deviation
// J. of the Amer. Statistical Assoc. (Theory and Methods), 88 (1993),p. 1273,
inline double Q_estimator(valarray<double>* err, size_t s)
{
	valarray<double> t1(binomial(s, 2));

	size_t i = 0;
	for (size_t w = 0; w < s - 1; w++)
	{
		for (size_t k = w + 1; k < s; k++)
		{
			t1[i] = (fabs((*err)[w] - (*err)[k]));
			i++;
		}
	}
	size_t h = s / 2 + 1;
	size_t k = binomial(h, 2) - 1;

#if __cplusplus == 201703L
	std::nth_element(std::execution::par, std::begin(t1), std::begin(t1) + k, std::end(t1));
#else
	std::nth_element(std::begin(t1), std::begin(t1) + k, std::end(t1));
#endif


	double cn[] = { 0,0,0.399, 0.994, 0.512 ,0.844 ,0.611, 0.857, 0.669 ,0.872 };
	double cc = 1;
	if (s <= 9)
	{
		cc = cn[s];
	}
	else
	{
		if (s % 2 != 0)
		{
			cc = s / (s + 1.4);
		}
		else
		{
			cc = s / (s + 3.8);
		}
	}
	return cc * 2.2219 * t1[k];
}

//computes the S estimator of an array of values with size s with a naive algorithm.
// The S estimator was published in
// Peter J. Rousseeuw, Christophe Croux, Alternatives to the Median-Absolute Deviation
// J. of the Amer. Statistical Assoc. (Theory and Methods), 88 (1993),p. 1273,

inline double S_estimator(valarray<double>* err, size_t s)
{
	valarray<double>m1(s);

	size_t sminus = s - 1;

	valarray<double>t1(sminus);

	for (size_t i = 0; i < s; i++)
	{
		size_t q = 0;
		for (size_t k = 0; k < s; k++)
		{
			if (i != k)
			{
				t1[q] = (fabs((*err)[i] - (*err)[k]));
				q++;
			}
		}
		m1[i] = (lowmedian(&t1, sminus));
	}
	double c = 1;
	double cn[] = { 0,0, 0.743, 1.851, 0.954 ,1.351, 0.993, 1.198 ,1.005, 1.131 };
	if (s <= 9)
	{
		c = cn[s];
	}
	else
	{
		if (s % 2 != 0)
		{
			c = s / (s - 0.9);
		}
	}

	return (c * 1.1926 * lowmedian(&m1, s));
}


//Computes the MAD estimator of an array and its median of values. The array should have a size size s, and a size/2 shalf 
inline double MAD_estimator(valarray<double>* err, double* m, size_t s, size_t shalf)
{
	*m = median(err, s, shalf);
	valarray<double>m1(s);
	for (size_t i = 0; i < s; i++)
	{
		m1[i] = (fabs((*err)[i] - *m));
	}
	double cn[] = { 0, 0, 1.196 ,1.495 ,1.363, 1.206, 1.200, 1.140, 1.129,1.107 };
	double c = 1.0;
	if (s <= 9)
	{
		c = cn[s];
	}
	else
	{
		c = s / (s - 0.8);
	}

	return  c * 1.4826 * median(&m1, s, shalf);
}

// Computes the biweightmidvariance of an array in one iteration. It also returns the median. It expects the size of the array and half the size of the array
// The biweight midvariance estimator  was described in 
// T. C. Beers,K. Flynn and K. Gebhardt,  Astron. J. 100 (1),32 (1990)
inline double onestepbiweightmidvariance(valarray<double>* err, double* median, size_t s, size_t shalf)
{
	double mad = MAD_estimator(err, median, s, shalf);
	double p1 = 0.0, p2 = 0.0;
	for (size_t i = 0; i < s; i++)
	{
		double ui = ((*err)[i] - *median) / (9.0 * mad);
		if (fabs(ui) < 1.0)
		{
			double ui2 = ui * ui;
			double g1 = (*err)[i] - *median;
			double g2 = (1.0 - ui2);
			p1 += g1 * g1 * pow(g2, 4.0);
			p2 += g2 * (1.0 - 5.0 * ui2);
		}
	}

	return (sqrt(s) * sqrt(p1) / fabs(p2));

}

inline double percentagebasedmidvariance(valarray<double>* err, double* med, size_t s, size_t shalf)
{
	double beta = 0.2;
	size_t m =(size_t) std::floor((1.0 - beta) * s + 0.5);
	*med = median(err,  s,  shalf);
	valarray<double> w(s);
	for (size_t i = 0; i < s; i++)
	{
		w[i] = fabs((*err)[i] - *med);
	}
	std::sort(std::begin(w), std::end(w));
	double wm = w[m];

	double sumpsi = 0,sumai=0;
	for (size_t i = 0; i < s; i++)
	{
		double Yi=((*err)[i] - *med) / wm;
		double h1 = (std::max(-1.0, std::min(1.0, Yi)));
		sumpsi +=h1*h1;
		sumai+= std::log((*err)[i]);
	}

	return (((double)s)*wm*sumpsi/(sumai*sumai));

}


//computes the T estimator for a array of size s. 
// The T estimator was described in Peter J. Rousseeuw, Christophe Croux, Alternatives to the Median-Absolute Deviation
// J. of the Amer. Statistical Assoc. (Theory and Methods), 88 (1993),p. 1273,
inline double T_estimator(valarray<double>* err, size_t s)
{

	valarray<double>med1(s);
	valarray<double>arr(s - 1);

	for (size_t i = 0; i < s; i++)
	{
		size_t q = 0;

		for (size_t j = 0; j < s; j++)
		{
			if (i != j)
			{
				arr[q] = (fabs((*err)[i] - (*err)[j]));
				q++;
			}
		}
		med1[i] = (median(&arr, s - 1, (s - 1) / 2));
	}
	size_t h = s / 2 + 1;
	double w = 0;
#if __cplusplus == 201703L
	sort(std::execution::par, std::begin(med1), std::end(med1));
#else
	sort(std::begin(med1), std::end(med1));
#endif

	for (size_t i = 1; i <= h; i++)
	{
		w += med1[i - 1];
	}
	return  (1.38 / ((double)h)) * w;
}


//expects a set of points with coordinates in (x,line_y), a boolean mask usedpoint which has values true for the indices in x and line_y to use for the fit. Then the function searches for the best hyperbolic fit
//with a minimum between minfocus and maxfocus, where the interval may be extended if the scale parameter is larger than 1. 
// the function then searches for points for which the usedpoint mask was set to false, that may be nevertheless included in the fit. It dpes this by comparing the errors of these points with supplied estimators and tolerance parameters
//the function then returns the modified mask.
inline bool findmodel(valarray<long>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, valarray<bool>* usedpoint, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method, size_t pointnumber, size_t pointnumberhalf)
{
	long this_focpos = 0;
	double thiserr = DBL_MAX, thisslope = 0.0, thisintercept = 0.0;
	//makes a fit with the initial usedpoint array as a mask
	if (!Search_min_error(x, line_y, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &this_focpos, use_median_regression, usedpoint))
	{
		return false;
	}

	vector<maybe_inliner> mp;
	mp.reserve(pointnumber);

	valarray<double>err(pointnumber);
	valarray<double>err_sqr(pointnumber);

	//makes a loop over all points
	for (size_t p = 0; p < pointnumber; p++)
	{

		double xh = ((double)(*x)[p] - (double)this_focpos);
		xh *= xh;
		//measures the error between the initial fit and point p
		double z = fabs(thisslope * xh + thisintercept) - (*line_y)[p];

		//store the error and the squared error
		err[p] = z;
		err_sqr[p] = z * z;

		//if point p is removed by the mask, add it as a maybe inlier
		if (!(*usedpoint)[p])
		{
			maybe_inliner o;
			o.point = p;
			o.error = z;
			mp.push_back(o);
		}
	}

	//if we have maybe inliers
	if (mp.size() > 0)
	{
		double m = 0, MAD = 0, average = 0, stdev = 0, S = 0, Q = 0, T = 0, biweightmidvariance = 0;

		//compute the user selected estimator for all the errors
		switch (rejection_method)
		{
		case tolerance_multiplies_standard_deviation_of_error:
		{
			stdeviation(&err, &stdev, &average, pointnumber);
			break;
		}
		case tolerance_is_significance_in_Grubbs_test:
		{
			stdeviation(&err, &stdev, &average, pointnumber);
			break;
		}
		case tolerance_is_decision_in_MAD_ESTIMATION:
		{
			MAD = MAD_estimator(&err, &m, pointnumber, pointnumberhalf);
			break;
		}
		case tolerance_is_biweight_midvariance:
		{
			biweightmidvariance = onestepbiweightmidvariance(&err, &m, pointnumber, pointnumberhalf);
			break;
		}
		case tolerance_is_percentagebased_midvariance:
		{
			biweightmidvariance = percentagebasedmidvariance(&err, &m, pointnumber, pointnumberhalf);
			break;
		}
		case tolerance_is_decision_in_S_ESTIMATION:
		{
			S = S_estimator(&err, pointnumber);
			m = median(&err, pointnumber, pointnumberhalf);
			break;
		}
		case tolerance_is_decision_in_Q_ESTIMATION:
		{
			m = median(&err, pointnumber, pointnumberhalf);
			Q = Q_estimator(&err, pointnumber);
			break;
		}
		case tolerance_is_decision_in_T_ESTIMATION:
		{
			m = median(&err, pointnumber, pointnumberhalf);
			T = T_estimator(&err, pointnumber);
			break;
		}
		case use_peirce_criterion:
		{
			stdeviation(&err_sqr, NULL, &average, pointnumber);
			break;
		}
		}

		//loop through the maybe inliers 
		for (size_t j = 0; j < mp.size(); j++)
		{
			//assume a point j is not an outlier
			bool isoutlier = false;

			// check if point j is within the bounds for an inlier with respect to the user selected estimator. If not, call it an outlier
			switch (rejection_method)
			{
			case tolerance_is_maximum_squared_error:
			{
				double G = mp[j].error * mp[j].error;
				if (G > fabs(tolerance))
				{
					isoutlier = true;
				}
				break;
			}
			case tolerance_multiplies_standard_deviation_of_error:
			{
				double G = fabs(mp[j].error - average);
				if (G > tolerance * stdev)
				{
					isoutlier = true;
				}
				break;
			}
			case tolerance_is_significance_in_Grubbs_test:
			{
				double G = fabs(mp[j].error - average);
				if (G > additionaldata * stdev)
				{
					isoutlier = true;
				}
				break;
			}
			case tolerance_is_decision_in_MAD_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / MAD);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case tolerance_is_biweight_midvariance:
			{
				double G = fabs((mp[j].error - m) / biweightmidvariance);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case tolerance_is_percentagebased_midvariance:
			{
				double G = fabs((mp[j].error - m) / biweightmidvariance);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}

			case tolerance_is_decision_in_S_ESTIMATION:
			{
				double G = (fabs(mp[j].error - m) / S);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}

			case tolerance_is_decision_in_Q_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / Q);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case tolerance_is_decision_in_T_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / T);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case use_peirce_criterion:
			{
				double G = mp[j].error * mp[j].error;
				if (G > average * additionaldata)
				{
					isoutlier = true;
				}
				break;
			}
			}
			//if point j is an inlier, modify the usedpoint bitmask
			if (isoutlier)
			{
				(*usedpoint)[mp[j].point] = false;
			}
			else
			{
				(*usedpoint)[mp[j].point] = true;
			}
		}
	}
	return true;
}
//the function which is exported by the library. It uses the ransac to compute a hyperbolic fit where it selects points as inliers and outliers. The documentation is provided in the header file.

bool focusposition_Regression(vector<long> x, vector<double> y, long* focpos, double* main_error, double* main_slope, double* main_intercept,
	vector<size_t>* indices_of_used_points,
	vector<double>* usedpoints_line_x, vector<double>* usedpoints_line_y, vector<size_t>* indices_of_removedpoints, vector<double>* removedpoints_line_x, vector<double>* removedpoints_line_y,
	double stop_after_seconds, size_t stop_after_numberofiterations_without_improvement, long backslash, double scale, bool use_median_regression,
	size_t maximum_number_of_outliers, outlier_criterion rejection_method, double tolerance)
{
	//test if garbage data was supplied.
	if (focpos == NULL)
	{
		return false;
	}
	if (x.size() < 4)
	{
		return false;
	}
	if (y.size() < x.size())
	{
		return false;
	}
	//set some often used variables
	size_t pointnumber = x.size();
	size_t pointnumberhalf = pointnumber / 2;
	size_t minimummodelsize = pointnumber - maximum_number_of_outliers;

	//make the supplied parameters consistent with each other in case there is a problem
	if (maximum_number_of_outliers == 0)
	{
		rejection_method = no_rejection;
	}

	if (rejection_method == no_rejection)
	{
		maximum_number_of_outliers = 0;
		minimummodelsize = pointnumber;
	}

	//ensure, we fit at least 4 points
	if (minimummodelsize < 4)
	{
		minimummodelsize = 4;
		maximum_number_of_outliers = pointnumber - minimummodelsize;
	}

	//test for a useless tolerance parameter with the Grubbs estimator
	if (rejection_method == tolerance_is_significance_in_Grubbs_test)
	{
		if (tolerance >= 0.9)
		{
			rejection_method = no_rejection;
		}
	}

	long minfocus = LONG_MAX, maxfocus = LONG_MIN;


	valarray<bool> indices(false,pointnumber);

	valarray<bool> indices2(pointnumber);
	valarray<double>line_yv(pointnumber);

	//set an initial bitmask, square the hfd values since we fit a hyperbola as a line where the hfd is squared, find the minimum and maximum focus value

	for (size_t i = 0; i < pointnumber; i++)
	{
		if (i >= maximum_number_of_outliers)
		{
			indices[i] = true;
		}

		line_yv[i] = y[i]*y[i];
		if (x[i] < minfocus)
		{
			minfocus = x[i];
		}
		if (x[i] > maxfocus)
		{
			maxfocus = x[i];
		}
	}

	//check if garbage data was supplied
	if (maxfocus - minfocus == 0)
		return false;


	double  error = DBL_MAX, intercept = 0, slope = 0, additionaldata = 0;

	//generate additional data for the Peirce criterion
	if (rejection_method == use_peirce_criterion)
	{
		additionaldata = peirce(pointnumber, maximum_number_of_outliers, 3);
	}
	//generate additional data for the Grubbs test
	if (rejection_method == tolerance_is_significance_in_Grubbs_test)
	{
		additionaldata = crit(tolerance, pointnumber);
	}


	valarray <long> xv(x.data(), x.size());

	//how many computations do we have to make based on the minimum fit model size
	size_t numbercomp = binomial(pointnumber, maximum_number_of_outliers);


#if __cplusplus == 201703L
	std::mutex mtx;
#endif


	//set some number of computations which should be distributed on several processors 
	const size_t number_of_attempts = 705432;


	//check if the supplied stop_after_numberofiterations_without_improvement was too small
	if (stop_after_numberofiterations_without_improvement < number_of_attempts)
	{
		stop_after_numberofiterations_without_improvement = number_of_attempts;
	}

	//if we have only a small number of computations to make
	if (numbercomp <= number_of_attempts)
	{
		std::unordered_set<std::vector<bool>> helper3;
		vector<valarray<bool>> arr(numbercomp);
		vector<valarray<bool>> arr2(numbercomp);
		//generate possible bitmasks
		for (size_t i = 0; i < numbercomp; i++)
		{
			arr[i] = indices;
			std::next_permutation(std::begin(indices), std::end(indices));
		}

		size_t count2 = 0;
		// the following code differs a bit in the syntax, if open-mp or C++17 is used. But the function is the same.
		// for each initial bitmask, call findmodel to fit the points to the bitmask and the bitmask if additional points need to be included.
		// copy the bitmask to a vector, which has a hash function and insert it to a list in which all bitmasks are inserted. If the insertion was successfull, the bitmask is new
		// it is then inserted into an array arr2 of different bitmasks.
		// increase the count of different bitmasks that were found
		// Let the initial error be at DBL_MAX.
		// then make a loop over the different final bitmasks we have found and make a curve fit for them. 
		// Store the parameters of the model if it has a smaller error than the model before.
#if __cplusplus == 201703L

		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
				if (b1)
				{
					std::vector<bool> helper(pointnumber);
					helper.assign(std::begin(arri), std::end(arri));
					bool b2;
					mtx.lock();
					b2 = helper3.insert(helper).second;
					if (b2)
					{
						arr2[count2] = arri;
						count2++;
					}
					mtx.unlock();
				}
			});
		std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
				if (b3)
				{
					mtx.lock();
					if (thiserr < error)
					{
						error = thiserr;
						*focpos = thisfocpos;
						slope = thisslope;
						intercept = thisintercept;
						indices2 = arri;
					}
					mtx.unlock();
				}
			});
#else
#pragma omp parallel for
		for (long i = 0; i < (long)numbercomp; i++)
		{
			bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
			if (b1)
			{
				std::vector<bool> helper(pointnumber);
				helper.assign(std::begin(arr[i]), std::end(arr[i]));
				bool b2;
#pragma omp critical
				{
					b2 = helper3.insert(helper).second;
					if (b2)
					{
						arr2[count2] = arr[i];
						count2++;
					}
				}
			}
		}

#pragma omp parallel for
		for (long i = 0; i < (long)count2; i++)
		{
			double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
			long thisfocpos = 0;
			bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
			if (b3)
			{
#pragma omp critical
				{
					if (thiserr < error)
					{
						error = thiserr;
						*focpos = thisfocpos;
						slope = thisslope;
						intercept = thisintercept;
						indices2 = arr2[i];
					}
				}
			}
		}
#endif
	}
	//if the number of possible bitmasks is just smaller than 100* number_of_attempts
	else if (numbercomp <= number_of_attempts * 100)
	{

		size_t p = (size_t)numbercomp / number_of_attempts;

		//make the same procedure as in the case with the smaller point number, but for p times. I.e. after the work was distributed to the processors start with another cycle and distribute again with
		// different bitmasks
		for (size_t o = 0; o < p; o++)
		{
			std::unordered_set<std::vector<bool>> helper3;
			vector<valarray<bool>> arr(number_of_attempts);
			vector<valarray<bool>> arr2(number_of_attempts);
			for (size_t i = 0; i < number_of_attempts; i++)
			{
				arr[i] = indices;
				std::next_permutation(std::begin(indices), std::end(indices));
			}

			// the following code differs a bit in the syntax, if open-mp or C++17 is used. But the function is the same.
			// for each initial bitmask, call findmodel to fit the points to the bitmask and the bitmask if additional points need to be included.
			// copy the bitmask to a vector, which has a hash function and insert it to a list in which all bitmasks are inserted. If the insertion was successfull, the bitmask is new
			// it is then inserted into an array arr2 of different bitmasks.
			// increase the count of different bitmasks that were found
			// Let the initial error be at DBL_MAX.
			// then make a loop over the different final bitmasks we have found and make a curve fit for them. 
			// Store the parameters of the model if it has a smaller error than the model before.

			size_t count2 = 0;
#if __cplusplus == 201703L
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
					if (b1)
					{
						std::vector<bool> helper(pointnumber);
						helper.assign(std::begin(arri), std::end(arri));
						bool b2;
						mtx.lock();
						b2 = helper3.insert(helper).second;
						if (b2)
						{
							arr2[count2] = arri;
							count2++;
						}
						mtx.unlock();
					}
				});

			std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
				{
					double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
					long thisfocpos = 0;
					bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
					if (b3)
					{
						mtx.lock();
						if (thiserr < error)
						{
							error = thiserr;
							*focpos = thisfocpos;
							slope = thisslope;
							intercept = thisintercept;
							indices2 = arri;
						}
						mtx.unlock();
					}
				});
#else			
#pragma omp parallel for
			for (long i = 0; i < (long)number_of_attempts; i++)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
				if (b1)
				{
					std::vector<bool> helper(pointnumber);
					helper.assign(std::begin(arr[i]), std::end(arr[i]));
					bool b2;
#pragma omp critical
					{
						b2 = helper3.insert(helper).second;
						if (b2)
						{
							arr2[count2] = arr[i];
							count2++;
						}
					}
				}
			}

#pragma omp parallel for
			for (long i = 0; i < count2; i++)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
				if (b3)
				{
#pragma omp critical
					{
						if (thiserr < error)
						{
							error = thiserr;
							*focpos = thisfocpos;
							slope = thisslope;
							intercept = thisintercept;
							indices2 = arr2[i];
						}
					}
				}
			}
#endif
		}

		std::unordered_set<std::vector<bool>> helper3;
		size_t s = (numbercomp % number_of_attempts) + 1;
		vector<valarray<bool>> arr(s);
		vector<valarray<bool>> arr2(s);

		//now distribute the work one last time for the remaining number of bitmasks.
		for (size_t i = 0; i < s; i++)
		{
			arr[i] = indices;
			std::next_permutation(std::begin(indices), std::end(indices));
		}




		size_t count2 = 0;

		// the following code differs a bit in the syntax, if open-mp or C++17 is used. But the function is the same.
		// for each initial bitmask, call findmodel to fit the points to the bitmask and the bitmask if additional points need to be included.
		// copy the bitmask to a vector, which has a hash function and insert it to a list in which all bitmasks are inserted. If the insertion was successfull, the bitmask is new
		// it is then inserted into an array arr2 of different bitmasks.
		// increase the count of different bitmasks that were found
		// Let the initial error be at DBL_MAX.
		// then make a loop over the different final bitmasks we have found and make a curve fit for them. 
		// Store the parameters of the model if it has a smaller error than the model before.
#if __cplusplus == 201703L
		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
				if (b1)
				{
					std::vector<bool> helper(pointnumber);
					helper.assign(std::begin(arri), std::end(arri));
					bool b2;
					mtx.lock();
					b2 = helper3.insert(helper).second;
					if (b2)
					{
						arr2[count2] = arri;
						count2++;
					}
					mtx.unlock();
				}
			});

		std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
				if (b3)
				{
					mtx.lock();
					if (thiserr < error)
					{
						error = thiserr;
						*focpos = thisfocpos;
						slope = thisslope;
						intercept = thisintercept;
						indices2 = arri;
					}
					mtx.unlock();
				}
			});
#else
#pragma omp parallel for
		for (long i = 0; i < (long)s; i++)
		{
			bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
			if (b1)
			{
				std::vector<bool> helper(pointnumber);
				helper.assign(std::begin(arr[i]), std::end(arr[i]));
				bool b2;
#pragma omp critical
				{
					b2 = helper3.insert(helper).second;
					if (b2)
					{
						arr2[count2] = arr[i];
						count2++;
					}
				}
			}
		}

#pragma omp parallel for
		for (long i = 0; i < count2; i++)
		{
			double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
			long thisfocpos = 0;
			bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
			if (b3)
			{
#pragma omp critical
				{
					if (thiserr < error)
					{
						error = thiserr;
						*focpos = thisfocpos;
						slope = thisslope;
						intercept = thisintercept;
						indices2 = arr2[i];
					}
				}
			}
		}
#endif
	}


	else
	{
		//we have too many points. So we start a true, non-deterministic ransac.
		std::random_device rng;
		std::mt19937 urng(rng());
		double seconds = 0;

#if __cplusplus == 201703L
		std::atomic<int> counter1{ 0 };
#else
		size_t counter1 = 0;
#endif

		auto start = std::chrono::steady_clock::now();

		do
		{
			//generate number_of_attempts random bitmasks and the helper structures
			std::unordered_set<std::vector<bool>> helper3;
			valarray<valarray<bool>> arr(number_of_attempts);
			valarray<valarray<bool>> arr2(number_of_attempts);
			for (size_t i = 0; i < number_of_attempts; i++)
			{
				shuffle(std::begin(indices), std::end(indices), urng);
				arr[i] = indices;
			}
			size_t count2 = 0;

			// the following code differs a bit in the syntax, if open-mp or C++17 is used. But the function is the same.
			// for each initial bitmask, call findmodel to fit the points to the bitmask and the bitmask if additional points need to be included.
			// copy the bitmask to a vector, which has a hash function and insert it to a list in which all bitmasks are inserted. If the insertion was successfull, the bitmask is new
			// it is then inserted into an array arr2 of different bitmasks.
			// increase the count of different bitmasks that were found
			// Let the initial error be at DBL_MAX.
			// then make a loop over the different final bitmasks we have found and make a curve fit for them. 
			// Store the parameters of the model if it has a smaller error than the model before.
#if __cplusplus == 201703L
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
					if (b1)
					{
						std::vector<bool> helper(pointnumber);
						helper.assign(std::begin(arri), std::end(arri));
						bool b2;
						mtx.lock();
						b2 = helper3.insert(helper).second;
						if (b2)
						{
							arr2[count2] = arri;
							count2++;
						}
						mtx.unlock();
					}
					//for each model we tried, increase counter1 by 1.
					counter1++;
				});

			std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
				{

					std::vector<bool> helper;
					helper.assign(std::begin(arri), std::end(arri));
					bool b2;
					mtx.lock();
					b2 = helper3.insert(helper).second;
					mtx.unlock();

					if (b2)
					{
						double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
						long thisfocpos = 0;
						bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
						if (b3)
						{
							mtx.lock();
							if (thiserr < error)
							{
								error = thiserr;
								*focpos = thisfocpos;
								slope = thisslope;
								intercept = thisintercept;
								//if we had found a better model set counter 1 =0
								counter1 = 0;
								indices2 = arri;
							}
							mtx.unlock();
						}
					}
				});
#else
#pragma omp parallel for
			for (long i = 0; i < number_of_attempts; i++)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
				//for each model we tried to find increase counter1 by 1
#pragma omp atomic
				counter1++;

				if (b1)
				{
					std::vector<bool> helper(pointnumber);
					helper.assign(std::begin(arr[i]), std::end(arr[i]));
					bool b2;
#pragma omp critical
					{
						b2 = helper3.insert(helper).second;
						if (b2)
						{
							arr2[count2] = arr[i];
							count2++;
						}
					}
				}
			}

#pragma omp parallel for
			for (long i = 0; i < count2; i++)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
				if (b3)
				{
#pragma omp critical
					{
						if (thiserr < error)
						{
							error = thiserr;
							*focpos = thisfocpos;
							slope = thisslope;
							intercept = thisintercept;
							indices2 = arr2[i];
							//if we had found a better model set the counter 1 to zero
							counter1 = 0;
						}
					}
				}
			}
#endif
			//if we had not found a better model after some user specified number of tries, stop the algorithm
			if (counter1 >= stop_after_numberofiterations_without_improvement)
			{
				break;
			}
			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			seconds = elapsed_seconds.count();
			//after the time specified by the user has passed, stop the algorithm
		} while (seconds < fabs(stop_after_seconds));
	}

	//if we had not found a good model return false
	if (error == DBL_MAX)
	{
		return false;
	}

	//correct the minimum focus for the user supplied backslash
	*focpos += backslash;

	//return the other parameters if specified
	if (main_error != NULL)
	{
		*main_error = error;
	}

	if (main_slope != NULL)
	{
		*main_slope = slope;
	}

	if (main_intercept != NULL)
	{
		*main_intercept = intercept;
	}

	//clear the arrays that were supplied as pointers to fill in return data.
	if (indices_of_removedpoints != NULL)
	{
		(*indices_of_removedpoints).clear();
	}

	if (indices_of_used_points != NULL)
	{
		(*indices_of_used_points).clear();
	}

	if (usedpoints_line_x != NULL)
	{
		(*usedpoints_line_x).clear();
	}

	if (usedpoints_line_y != NULL)
	{
		(*usedpoints_line_y).clear();
	}

	if (removedpoints_line_x != NULL)
	{
		(*removedpoints_line_x).clear();
	}

	if (removedpoints_line_y != NULL)
	{
		(*removedpoints_line_y).clear();
	}

	//fill in the return data for the indices of the used points and the coordinates of the fitted model in various formats (as a line and as a hyperbola).
	for (size_t k = 0; k < pointnumber; k++)
	{
		if (indices2[k])
		{
			if (indices_of_used_points != NULL)
			{
				(*indices_of_used_points).push_back(k);
			}

			if (usedpoints_line_x != NULL)
			{
				double k1 = (double)x[k] - (double)*focpos;
				k1 *= k1;
				(*usedpoints_line_x).push_back(k1);
			}

			if (usedpoints_line_y != NULL)
			{
				double w = y[k];
				(*usedpoints_line_y).push_back(w * w);
			}
		}
		else
		{
			if (indices_of_removedpoints != NULL)
			{
				(*indices_of_removedpoints).push_back(k);
			}

			if (removedpoints_line_x != NULL)
			{
				double k1 = (double)x[k] - (double)*focpos;
				k1 *= k1;
				(*removedpoints_line_x).push_back(k1);
			}

			if (removedpoints_line_y != NULL)
			{
				double w = y[k];
				(*removedpoints_line_y).push_back(w * w);
			}
		}

	}


	return true;
}


//the function is exported by the library. It starts 2 hyperbolic fits from 2 point sets and computes the difference between the minima of the hyperbolas. 
//The documentation of this function is provided in the header file.

bool findbackslash_Regression(long* backslash,
	vector<long> x1, vector<double> y1, vector<long> x2, vector<double> y2, double* main_error1, double* main_slope1, double* main_intercept1, vector<size_t>* indicesofusedpoints1,
	vector<double>* used_points1_line_x, vector<double>* used_points1_line_y, vector<size_t>* indicesofremovedpoints1, vector<double>* removedpoints1_line_x, vector<double>* removedpoints1_line_y,
	double* main_error2, double* main_slope2, double* main_intercept2, vector<size_t>* indicesofusedpoints2,
	vector<double>* used_points2_line_x, vector<double>* used_points2_line_y, vector<size_t>* indicesofremovedpoints2, vector<double>* removedpoints2_line_x, vector<double>* removedpoints2_line_y,
	double stop_after_seconds, size_t stop_after_numberofiterations_without_improvement, double scale, bool use_median_regression,
	size_t maximum_number_of_outliers, outlier_criterion rejection_method, double tolerance)
{

	if (backslash == NULL)
	{
		return false;
	}

	long focpos1 = 0, focpos2 = 0;
	//make the first fit;
	if (!focusposition_Regression(x1, y1, &focpos1, main_error1, main_slope1, main_intercept1, indicesofusedpoints1, used_points1_line_x, used_points1_line_y, indicesofremovedpoints1, removedpoints1_line_x, removedpoints1_line_y, stop_after_seconds,
		stop_after_numberofiterations_without_improvement, 0, scale, use_median_regression, maximum_number_of_outliers, rejection_method, tolerance))
	{
		return false;
	}
	//make the second fit;
	if (!focusposition_Regression(x2, y2, &focpos2, main_error2, main_slope2, main_intercept2, indicesofusedpoints2, used_points2_line_x, used_points2_line_y, indicesofremovedpoints2, removedpoints2_line_x, removedpoints2_line_y, stop_after_seconds,
		stop_after_numberofiterations_without_improvement, 0, scale, use_median_regression, maximum_number_of_outliers, rejection_method, tolerance))
	{
		return false;
	}
	//return the computed backslash
	*backslash = focpos2 - focpos1;
	return true;
}

//returns the pstatus variable. status=0 means the image class was successfully constructed.
int image::status()
{
	return pstatus;
}
//returns the power of an image in fourier mode.
double image::invpower()
{
	return pinvpower;
}


//returns the focuser position of the image class
long image::focuser_position()
{
	return pfocuser_position;
}

//returns the hfd that can be supplied in the constructor(note that this value is not used in the algorithms, but can be stored as an option, for applications which
// have their own hdf analysis algorithm.)
double image::hfd()
{
	return phfd;
}
//returns the fwhm (note that this value is not used in the algorithms, but can be stored as an option, for applications which
// have their own hdf analysis algorithm.)
double image::fwhm()
{
	return pfwhm;
}

//constructs the image class from a fits file given by filename. The focuser position of the image is either supplied as a
//parameter or read from the fits file. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(string* filename, long focuser_position, double hfd, double fwhm)
{
	pstatus = 0;
	phfd = hfd;
	pfwhm = fwhm;
	pfocuser_position = focuser_position;
	fitsfile* fptr;
	fits_open_file(&fptr, filename[0].c_str(), READONLY, &pstatus);
	if (pstatus !=0)
	{
		throw std::invalid_argument("could not open the fits file");
	}
	else
	{
		pinvpower=datafilling0(fptr, focuser_position);
		fits_close_file(fptr, &pstatus);
		if (pstatus != 0)
		{
			throw std::invalid_argument("could not close the fits file");
		}
	}
}
// constructs the image class from a fits file object. The focuser position of the image is either supplied as a
// parameter or read from the fits file. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(fitsfile * fptr, long focuser_position, double hfd , double fwhm )
{
	pstatus = 0;
	phfd = hfd;
	pfwhm = fwhm;
	pfocuser_position = focuser_position;
	pinvpower = datafilling0(fptr, focuser_position);
}

//fills the focuser position and power function, and optionally hfd and fwhm values into the class data.
//The focuser position may be supplied explicitely or given in the fits files under Keywords
//"FOCUSPOS", "FOCUSERPOS", "FOCUSERPOSITION", "FOCUSPOSITION".
double image::datafilling0(fitsfile* fptr, long focuser_position)
{

	size_t dimension1 = 0;
	size_t dimension2 = 0;
	fits_read_key(fptr, TLONG, "NAXIS1", &dimension1, NULL, &pstatus);
	if (pstatus != 0)
	{
		throw std::invalid_argument("could not read fits file dimension 1");
	}
	fits_read_key(fptr, TLONG, "NAXIS2", &dimension2, NULL, &pstatus);
	if (pstatus != 0)
	{
		throw std::invalid_argument("could not read fits file dimension 2");
	}
	if (dimension1 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension1 was <=2");
	}
	if (dimension2 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension2 was <=2");
	}
	int status = 0;
	if (focuser_position == LONG_MIN)
	{
		if (fits_read_key(fptr, TLONG, "FOCUSPOS", &pfocuser_position, NULL, &status) == KEY_NO_EXIST)
		{
			if (fits_read_key(fptr, TLONG, "FOCUSERPOS", &pfocuser_position, NULL, &status) == KEY_NO_EXIST)
			{
				if (fits_read_key(fptr, TLONG, "FOCUSERPOSITION", &pfocuser_position, NULL, &status) == KEY_NO_EXIST)
				{
					if (fits_read_key(fptr, TLONG, "FOCUSPOSITION", &pfocuser_position, NULL, &status) == KEY_NO_EXIST)
					{
						if (fits_read_key(fptr, TLONG, "FOCUSMOTORPOSITION", &pfocuser_position, NULL, &status) == KEY_NO_EXIST)
						{
							if (fits_read_key(fptr, TLONG, "FOCUSMOTORPOS", &pfocuser_position, NULL, &status) == KEY_NO_EXIST)
							{
								status = -1;
							}
							else
							{
								status = 0;
							}
						}
					}
				}
			}
		}
		if (status != 0)
		{
			pstatus = -1;
			throw std::invalid_argument("could not read the focuser position from the fits file");
		}
	}
	else
	{
		pfocuser_position = focuser_position;
		pstatus = 0;
	}
	int bitpix;
	fits_get_img_type(fptr, &bitpix, &pstatus);
	if (pstatus != 0)
	{
		throw std::invalid_argument("could not read fits file imae type");
	}
	if (bitpix == -64)
	{
		double nulval = NULL;
		int anynull;

		std::vector<double> pimagedata_double;
		pimagedata_double.resize(dimension1 * dimension2);

		fits_read_img(fptr, TDOUBLE, 1, dimension1 * dimension2, &nulval, &pimagedata_double[0], &anynull, &pstatus);
		if (pstatus != 0)
		{
			throw std::invalid_argument("could not read fits file image data");
		}
		else
		{
			return fouriertransform(NULL, &pimagedata_double, dimension1, dimension2);
		}
	}
	else
	{
		float nulval = NULL;
		int anynull;

		std::vector<float> pimagedata_float;
		pimagedata_float.resize(dimension1 * dimension2);

		fits_read_img(fptr, TFLOAT, 1, dimension1 * dimension2, &nulval, &pimagedata_float[0], &anynull, &pstatus);
		if (pstatus != 0)
		{
			throw std::invalid_argument("could not read fits file image data");
		}
		else
		{
			return fouriertransform(&pimagedata_float, NULL, dimension1, dimension2);
		}
	}

}

//analyzes the fourier data of the image and calculates the power. The image data can be either in fload or double.
double image::fouriertransform(vector<float>* p4, vector<double>* p5, size_t dimension1, size_t dimension2)
{

	if (p4 != NULL)
	{
		//compute the dft of the image
		Mat I(dimension2, dimension1, CV_32F, (*p4).data());

		Mat planes[] = { Mat_<float>(I), Mat::zeros(I.size(), 5) };
		Mat complexI;
		merge(planes, 2, complexI);

		dft(complexI, complexI);
		split(complexI, planes);
		// compute the magnitude
		Mat magI = planes[0].mul(planes[0]) + planes[1].mul(planes[1]);

		//sum over rows
		Mat sumcol;
		cv::reduce(magI, sumcol, 0, cv::REDUCE_SUM, -1);
		//compute the offset
		size_t avestart = (size_t)(round((magI.cols) / 2 - (magI.cols) / 8)) - 1;
		size_t length1 = (size_t)((magI.cols) / 2)-1 - avestart;
		Scalar offset = cv::mean(sumcol(cv::Rect(avestart, 0, length1, 1)));
		//subtract the offset and clip
		sumcol -= offset;
		sumcol = cv::max(sumcol, 0.0);
		//compute the power
		size_t start = 1;
		size_t length2 = (size_t)(magI.cols / 2) - start;

		Scalar power = sum(sumcol(cv::Rect(start, 0, length2, 1)));

		return 1.0/(power[0] + power[1] + power[2] + power[3]);
	}
	else
	{
		if (p5 == NULL)
		{
			pstatus = -1;
			throw std::invalid_argument("no valid input data");
			return -1;
		}
		//compute the fourier transform for a double image
		Mat I(dimension2, dimension1, CV_64F, (*p5).data());

		Mat planes[] = { Mat_<double>(I), Mat::zeros(I.size(), 6) };
		Mat complexI;
		merge(planes, 2, complexI);
		dft(complexI, complexI);
		split(complexI, planes);
		//compute the magnitude
		Mat magI = planes[0].mul(planes[0]) + planes[1].mul(planes[1]);
		//sum over rows
		Mat sumcol;
		cv::reduce(magI, sumcol, 0, cv::REDUCE_SUM, -1);
		//compute the offset
		size_t avestart = (size_t)(round((magI.cols) / 2 - (magI.cols) / 8)) - 1;
		size_t length1 = (size_t)((magI.cols) / 2)-1 - avestart;
		Scalar offset = cv::mean(sumcol(cv::Rect(avestart, 0, length1, 1)));
		//subtract the offset and clip
		sumcol -= offset;
		sumcol = cv::max(sumcol, 0.0);
		//compute the power
		size_t start = 1;
		size_t length2 = (size_t)(magI.cols / 2) - start;

		Scalar power = sum(sumcol(cv::Rect(start, 0, length2, 1)));

		return 1.0/( power[0] + power[1] + power[2] + power[3]);
	}

}

//constructs the image class from a double array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t dimension1, size_t dimension2, long focuser_position, vector<double> * imagedata, double hfd, double fwhm)
{
	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (dimension1 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension1 was <=2");
	}
	if (dimension2 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension2 was <=2");
	}
	if (imagedata->size() ==0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((dimension1 * dimension2) != imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}
	pinvpower=fouriertransform(NULL, imagedata, dimension1, dimension2);
}
//constructs the image class from a float array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t dimension1, size_t dimension2, long focuser_position, vector<float> * imagedata, double hfd, double fwhm)
{
	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (dimension1 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension1 was <=2");
	}
	if (dimension2 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension2 was <=2");
	}
	if (imagedata->size() == 0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((dimension1*dimension2)!= imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}
	pinvpower = fouriertransform(imagedata, NULL, dimension1, dimension2);
}
//constructs the image class from a long long array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t dimension1, size_t dimension2, long focuser_position, vector<long long> * imagedata, double hfd, double fwhm)
{
	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (dimension1 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension1 was <=2");
	}
	if (dimension2 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension2 was <=2");
	}
	if (imagedata->size() == 0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((dimension1 * dimension2) != imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}

	vector<float>floatvec(imagedata->begin(), imagedata->end());
	pinvpower = fouriertransform(&floatvec, NULL, dimension1, dimension2);
}
//constructs the image class from a long array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t dimension1, size_t dimension2, long focuser_position, vector<long> * imagedata, double hfd, double fwhm)
{

	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (dimension1 <=2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension1 was <=2");
	}
	if (dimension2 <=2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension2 was <=2");
	}
	if (imagedata->size() == 0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((dimension1 * dimension2) != imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}


	vector<float>floatvec(imagedata->begin(), imagedata->end());
	pinvpower = fouriertransform(&floatvec, NULL, dimension1, dimension2);
}
//constructs the image class from a short array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t dimension1, size_t dimension2, long focuser_position, vector<short> * imagedata, double hfd, double fwhm)
{
	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (dimension1 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension1 was <=2");
	}
	if (dimension2 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension2 was <=2");
	}
	if (imagedata->size() == 0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((dimension1 * dimension2) != imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}

	
	vector<float>floatvec(imagedata->begin(), imagedata->end());
	pinvpower = fouriertransform(&floatvec, NULL, dimension1, dimension2);
}
//constructs the image class from an int8_t array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t dimension1, size_t dimension2, long focuser_position, vector<int8_t> * imagedata, double hfd, double fwhm)
{

	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (dimension1 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension1 was <=2");
	}
	if (dimension2 <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("dimension2 was <=2");
	}
	if (imagedata->size() == 0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((dimension1 * dimension2) != imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}


	vector<float>floatvec(imagedata->begin(), imagedata->end());
	pinvpower = fouriertransform(&floatvec, NULL, dimension1, dimension2);
}

//computes the focus point from a vector of validly constructed image classes which must all have status 0.
//the other parameters are as in focusposition_Regression with 2 exceptions: double* main_slope, double* main_intercept are the slope 
// and intercept for the functions invpower=slope(x-focus_point)^2+intercept, where invpower is the given by the invpower method 
//if x is the focus motor position of an image class. It is the inverse of the power function from the fourier analysis.
//also, the parameters  from focusposition_Regression where the fitted curve is returned in a coordinate system where it is represented 
//by a line are omitted. 

bool focusposition_Regression2(std::vector<image>* images, long* focpos, double* main_error, double* main_slope, double* main_intercept, double* theta,
	vector<size_t>* indices_of_used_points, vector<size_t>* indices_of_removedpoints, double stop_after_seconds, size_t stop_after_numberofiterations_without_improvement, long backslash, double scale, bool use_median_regression,
	size_t maximum_number_of_outliers, outlier_criterion rejection_method, double tolerance)
{
	//test if garbage data was supplied.
	if (focpos == NULL)
	{
		return false;
	}
	size_t pointnumber = 0;
	for (size_t i = 0; i < (*images).size(); i++)
	{
		if ((*images)[i].status() == 0)
		{
			pointnumber++;
		}
	}
	if (pointnumber < 4)
	{
		return false;
	}

	//set some often used variables
	size_t pointnumberhalf = pointnumber / 2;
	size_t minimummodelsize = pointnumber - maximum_number_of_outliers;

	//make the supplied parameters consistent with each other in case there is a problem
	if (maximum_number_of_outliers == 0)
	{
		rejection_method = no_rejection;
	}

	if (rejection_method == no_rejection)
	{
		maximum_number_of_outliers = 0;
		minimummodelsize = pointnumber;
	}

	//ensure, we fit at least 4 points
	if (minimummodelsize < 4)
	{
		minimummodelsize = 4;
		maximum_number_of_outliers = pointnumber - minimummodelsize;
	}

	//test for a useless tolerance parameter with the Grubbs estimator
	if (rejection_method == tolerance_is_significance_in_Grubbs_test)
	{
		if (tolerance >= 0.9)
		{
			rejection_method = no_rejection;
		}
	}

	long minfocus = LONG_MAX, maxfocus = LONG_MIN;


	valarray<bool> indices(false, pointnumber);

	valarray<bool> indices2(pointnumber);
	valarray<double>line_yv(pointnumber);
	valarray <long> xv(pointnumber);
	//for all image classes which were successfully created, set an initial bitmask, set the y values to the inverse power, find the minimum and maximum focus value

	for (size_t i = 0; i < pointnumber; i++)
	{
		if ((*images)[i].status() == 0)
		{

			if (i >= maximum_number_of_outliers)
			{
				indices[i] = true;
			}
			xv[i] = (*images)[i].focuser_position();

			line_yv[i] =  (*images)[i].invpower();
			if (xv[i] < minfocus)
			{
				minfocus = xv[i];
			}
			if (xv[i] > maxfocus)
			{
				maxfocus = xv[i];
			}
		}
	}

	//check if garbage data was supplied
	if (maxfocus - minfocus == 0)
		return false;


	double  error = DBL_MAX, intercept = 0, slope = 0, additionaldata = 0;

	//generate additional data for the Peirce criterion
	if (rejection_method == use_peirce_criterion)
	{
		additionaldata = peirce(pointnumber, maximum_number_of_outliers, 3);
	}
	//generate additional data for the Grubbs test
	if (rejection_method == tolerance_is_significance_in_Grubbs_test)
	{
		additionaldata = crit(tolerance, pointnumber);
	}



	//how many computations do we have to make based on the minimum fit model size
	size_t numbercomp = binomial(pointnumber, maximum_number_of_outliers);


#if __cplusplus == 201703L
	std::mutex mtx;
#endif


	//set some number of computations which should be distributed on several processors 
	const size_t number_of_attempts = 705432;


	//check if the supplied stop_after_numberofiterations_without_improvement was too small
	if (stop_after_numberofiterations_without_improvement < number_of_attempts)
	{
		stop_after_numberofiterations_without_improvement = number_of_attempts;
	}

	//if we have only a small number of computations to make
	if (numbercomp <= number_of_attempts)
	{
		std::unordered_set<std::vector<bool>> helper3;
		vector<valarray<bool>> arr(numbercomp);
		vector<valarray<bool>> arr2(numbercomp);
		//generate possible bitmasks
		for (size_t i = 0; i < numbercomp; i++)
		{
			arr[i] = indices;
			std::next_permutation(std::begin(indices), std::end(indices));
		}

		size_t count2 = 0;
		// the following code differs a bit in the syntax, if open-mp or C++17 is used. But the function is the same.
		// for each initial bitmask, call findmodel to fit the points to the bitmask and the bitmask if additional points need to be included.
		// copy the bitmask to a vector, which has a hash function and insert it to a list in which all bitmasks are inserted. If the insertion was successfull, the bitmask is new
		// it is then inserted into an array arr2 of different bitmasks.
		// increase the count of different bitmasks that were found
		// Let the initial error be at DBL_MAX.
		// then make a loop over the different final bitmasks we have found and make a curve fit for them. 
		// Store the parameters of the model if it has a smaller error than the model before.
#if __cplusplus == 201703L

		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
				if (b1)
				{
					std::vector<bool> helper(pointnumber);
					helper.assign(std::begin(arri), std::end(arri));
					bool b2;
					mtx.lock();
					b2 = helper3.insert(helper).second;
					if (b2)
					{
						arr2[count2] = arri;
						count2++;
					}
					mtx.unlock();
				}
			});
		std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
				if (b3)
				{
					mtx.lock();
					if (thiserr < error)
					{
						error = thiserr;
						*focpos = thisfocpos;
						slope = thisslope;
						intercept = thisintercept;
						indices2 = arri;
					}
					mtx.unlock();
				}
			});
#else
#pragma omp parallel for
		for (long i = 0; i < (long)numbercomp; i++)
		{
			bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
			if (b1)
			{
				std::vector<bool> helper(pointnumber);
				helper.assign(std::begin(arr[i]), std::end(arr[i]));
				bool b2;
#pragma omp critical
				{
					b2 = helper3.insert(helper).second;
					if (b2)
					{
						arr2[count2] = arr[i];
						count2++;
					}
				}
			}
		}

#pragma omp parallel for
		for (long i = 0; i < (long)count2; i++)
		{
			double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
			long thisfocpos = 0;
			bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
			if (b3)
			{
#pragma omp critical
				{
					if (thiserr < error)
					{
						error = thiserr;
						*focpos = thisfocpos;
						slope = thisslope;
						intercept = thisintercept;
						indices2 = arr2[i];
					}
				}
			}
		}
#endif
	}
	//if the number of possible bitmasks is just smaller than 100* number_of_attempts
	else if (numbercomp <= number_of_attempts * 100)
	{

		size_t p = (size_t)numbercomp / number_of_attempts;

		//make the same procedure as in the case with the smaller point number, but for p times. I.e. after the work was distributed to the processors start with another cycle and distribute again with
		// different bitmasks
		for (size_t o = 0; o < p; o++)
		{
			std::unordered_set<std::vector<bool>> helper3;
			vector<valarray<bool>> arr(number_of_attempts);
			vector<valarray<bool>> arr2(number_of_attempts);
			for (size_t i = 0; i < number_of_attempts; i++)
			{
				arr[i] = indices;
				std::next_permutation(std::begin(indices), std::end(indices));
			}

			// the following code differs a bit in the syntax, if open-mp or C++17 is used. But the function is the same.
			// for each initial bitmask, call findmodel to fit the points to the bitmask and the bitmask if additional points need to be included.
			// copy the bitmask to a vector, which has a hash function and insert it to a list in which all bitmasks are inserted. If the insertion was successfull, the bitmask is new
			// it is then inserted into an array arr2 of different bitmasks.
			// increase the count of different bitmasks that were found
			// Let the initial error be at DBL_MAX.
			// then make a loop over the different final bitmasks we have found and make a curve fit for them. 
			// Store the parameters of the model if it has a smaller error than the model before.

			size_t count2 = 0;
#if __cplusplus == 201703L
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
					if (b1)
					{
						std::vector<bool> helper(pointnumber);
						helper.assign(std::begin(arri), std::end(arri));
						bool b2;
						mtx.lock();
						b2 = helper3.insert(helper).second;
						if (b2)
						{
							arr2[count2] = arri;
							count2++;
						}
						mtx.unlock();
					}
				});

			std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
				{
					double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
					long thisfocpos = 0;
					bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
					if (b3)
					{
						mtx.lock();
						if (thiserr < error)
						{
							error = thiserr;
							*focpos = thisfocpos;
							slope = thisslope;
							intercept = thisintercept;
							indices2 = arri;
						}
						mtx.unlock();
					}
				});
#else			
#pragma omp parallel for
			for (long i = 0; i < (long)number_of_attempts; i++)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
				if (b1)
				{
					std::vector<bool> helper(pointnumber);
					helper.assign(std::begin(arr[i]), std::end(arr[i]));
					bool b2;
#pragma omp critical
					{
						b2 = helper3.insert(helper).second;
						if (b2)
						{
							arr2[count2] = arr[i];
							count2++;
						}
					}
				}
			}

#pragma omp parallel for
			for (long i = 0; i < count2; i++)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
				if (b3)
				{
#pragma omp critical
					{
						if (thiserr < error)
						{
							error = thiserr;
							*focpos = thisfocpos;
							slope = thisslope;
							intercept = thisintercept;
							indices2 = arr2[i];
						}
					}
				}
			}
#endif
		}

		std::unordered_set<std::vector<bool>> helper3;
		size_t s = (numbercomp % number_of_attempts) + 1;
		vector<valarray<bool>> arr(s);
		vector<valarray<bool>> arr2(s);

		//now distribute the work one last time for the remaining number of bitmasks.
		for (size_t i = 0; i < s; i++)
		{
			arr[i] = indices;
			std::next_permutation(std::begin(indices), std::end(indices));
		}




		size_t count2 = 0;

		// the following code differs a bit in the syntax, if open-mp or C++17 is used. But the function is the same.
		// for each initial bitmask, call findmodel to fit the points to the bitmask and the bitmask if additional points need to be included.
		// copy the bitmask to a vector, which has a hash function and insert it to a list in which all bitmasks are inserted. If the insertion was successfull, the bitmask is new
		// it is then inserted into an array arr2 of different bitmasks.
		// increase the count of different bitmasks that were found
		// Let the initial error be at DBL_MAX.
		// then make a loop over the different final bitmasks we have found and make a curve fit for them. 
		// Store the parameters of the model if it has a smaller error than the model before.
#if __cplusplus == 201703L
		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
				if (b1)
				{
					std::vector<bool> helper(pointnumber);
					helper.assign(std::begin(arri), std::end(arri));
					bool b2;
					mtx.lock();
					b2 = helper3.insert(helper).second;
					if (b2)
					{
						arr2[count2] = arri;
						count2++;
					}
					mtx.unlock();
				}
			});

		std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
				if (b3)
				{
					mtx.lock();
					if (thiserr < error)
					{
						error = thiserr;
						*focpos = thisfocpos;
						slope = thisslope;
						intercept = thisintercept;
						indices2 = arri;
					}
					mtx.unlock();
				}
			});
#else
#pragma omp parallel for
		for (long i = 0; i < (long)s; i++)
		{
			bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
			if (b1)
			{
				std::vector<bool> helper(pointnumber);
				helper.assign(std::begin(arr[i]), std::end(arr[i]));
				bool b2;
#pragma omp critical
				{
					b2 = helper3.insert(helper).second;
					if (b2)
					{
						arr2[count2] = arr[i];
						count2++;
					}
				}
			}
		}

#pragma omp parallel for
		for (long i = 0; i < count2; i++)
		{
			double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
			long thisfocpos = 0;
			bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
			if (b3)
			{
#pragma omp critical
				{
					if (thiserr < error)
					{
						error = thiserr;
						*focpos = thisfocpos;
						slope = thisslope;
						intercept = thisintercept;
						indices2 = arr2[i];
					}
				}
			}
		}
#endif
	}


	else
	{
		//we have too many points. So we start a true, non-deterministic ransac.
		std::random_device rng;
		std::mt19937 urng(rng());
		double seconds = 0;

#if __cplusplus == 201703L
		std::atomic<int> counter1{ 0 };
#else
		size_t counter1 = 0;
#endif

		auto start = std::chrono::steady_clock::now();

		do
		{
			//generate number_of_attempts random bitmasks and the helper structures
			std::unordered_set<std::vector<bool>> helper3;
			valarray<valarray<bool>> arr(number_of_attempts);
			valarray<valarray<bool>> arr2(number_of_attempts);
			for (size_t i = 0; i < number_of_attempts; i++)
			{
				shuffle(std::begin(indices), std::end(indices), urng);
				arr[i] = indices;
			}
			size_t count2 = 0;

			// the following code differs a bit in the syntax, if open-mp or C++17 is used. But the function is the same.
			// for each initial bitmask, call findmodel to fit the points to the bitmask and the bitmask if additional points need to be included.
			// copy the bitmask to a vector, which has a hash function and insert it to a list in which all bitmasks are inserted. If the insertion was successfull, the bitmask is new
			// it is then inserted into an array arr2 of different bitmasks.
			// increase the count of different bitmasks that were found
			// Let the initial error be at DBL_MAX.
			// then make a loop over the different final bitmasks we have found and make a curve fit for them. 
			// Store the parameters of the model if it has a smaller error than the model before.
#if __cplusplus == 201703L
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
					if (b1)
					{
						std::vector<bool> helper(pointnumber);
						helper.assign(std::begin(arri), std::end(arri));
						bool b2;
						mtx.lock();
						b2 = helper3.insert(helper).second;
						if (b2)
						{
							arr2[count2] = arri;
							count2++;
						}
						mtx.unlock();
					}
					//for each model we tried, increase counter1 by 1.
					counter1++;
				});

			std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
				{

					std::vector<bool> helper;
					helper.assign(std::begin(arri), std::end(arri));
					bool b2;
					mtx.lock();
					b2 = helper3.insert(helper).second;
					mtx.unlock();

					if (b2)
					{
						double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
						long thisfocpos = 0;
						bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
						if (b3)
						{
							mtx.lock();
							if (thiserr < error)
							{
								error = thiserr;
								*focpos = thisfocpos;
								slope = thisslope;
								intercept = thisintercept;
								//if we had found a better model set counter 1 =0
								counter1 = 0;
								indices2 = arri;
							}
							mtx.unlock();
						}
					}
				});
#else
#pragma omp parallel for
			for (long i = 0; i < number_of_attempts; i++)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
				//for each model we tried to find increase counter1 by 1
#pragma omp atomic
				counter1++;

				if (b1)
				{
					std::vector<bool> helper(pointnumber);
					helper.assign(std::begin(arr[i]), std::end(arr[i]));
					bool b2;
#pragma omp critical
					{
						b2 = helper3.insert(helper).second;
						if (b2)
						{
							arr2[count2] = arr[i];
							count2++;
						}
					}
				}
			}

#pragma omp parallel for
			for (long i = 0; i < count2; i++)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
				if (b3)
				{
#pragma omp critical
					{
						if (thiserr < error)
						{
							error = thiserr;
							*focpos = thisfocpos;
							slope = thisslope;
							intercept = thisintercept;
							indices2 = arr2[i];
							//if we had found a better model set the counter 1 to zero
							counter1 = 0;
						}
					}
				}
			}
#endif
			//if we had not found a better model after some user specified number of tries, stop the algorithm
			if (counter1 >= stop_after_numberofiterations_without_improvement)
			{
				break;
			}
			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			seconds = elapsed_seconds.count();
			//after the time specified by the user has passed, stop the algorithm
		} while (seconds < fabs(stop_after_seconds));
	}

	//if we had not found a good model return false
	if (error == DBL_MAX)
	{
		return false;
	}

	//correct the minimum focus for the user supplied backslash
	*focpos += backslash;

	//return the other parameters if specified
	if (main_error != NULL)
	{
		*main_error = error;
	}

	if (main_slope != NULL)
	{
		*main_slope = slope;
	}

	if (main_intercept != NULL)
	{
		*main_intercept = intercept;
	}

	//clear the arrays that were supplied as pointers to fill in return data.
	if (indices_of_removedpoints != NULL)
	{
		(*indices_of_removedpoints).clear();
	}

	if (indices_of_used_points != NULL)
	{
		(*indices_of_used_points).clear();
	}



	//fill in the return data for the indices of the used points and the coordinates of the fitted model in various formats (as a line and as a hyperbola).
	for (size_t k = 0; k < pointnumber; k++)
	{
		if (indices2[k])
		{
			if (indices_of_used_points != NULL)
			{
				(*indices_of_used_points).push_back(k);
			}
		}
		else
		{
			if (indices_of_removedpoints != NULL)
			{
				(*indices_of_removedpoints).push_back(k);
			}
		}

	}


	return true;
}

bool findbackslash_Regression2(long* backslash,
	vector<image>*images1, vector<image>*images2, double* main_error1 , double* main_slope1 , double* main_intercept1 , vector<size_t>*indicesofusedpoints1, vector<size_t>*indicesofremovedpoints1 ,
	double* main_error2 , double* main_slope2 , double* main_intercept2 , double* theta1 , double* theta2 , vector<size_t>*indicesofusedpoints2 , vector<size_t>*indicesofremovedpoints2 ,
	double stop_after_seconds , size_t stop_after_numberofiterations_without_improvement , double scale , bool use_median_regression ,
	size_t maximum_number_of_outliers , outlier_criterion rejection_method , double tolerance)
{

	if (backslash == NULL)
	{
		return false;
	}

	long focpos1 = 0, focpos2 = 0;
	//make the first fit;
	if (!focusposition_Regression2(images1, &focpos1, main_error1, main_slope1, main_intercept1,theta1, indicesofusedpoints1, indicesofremovedpoints1,stop_after_seconds,
		stop_after_numberofiterations_without_improvement, 0, scale, use_median_regression, maximum_number_of_outliers, rejection_method, tolerance))
	{
		return false;
	}
	//make the second fit;
	if (!focusposition_Regression2(images2, &focpos2, main_error2, main_slope2, main_intercept2, theta2 ,indicesofusedpoints2, indicesofremovedpoints2, stop_after_seconds,
		stop_after_numberofiterations_without_improvement, 0, scale, use_median_regression, maximum_number_of_outliers, rejection_method, tolerance))
	{
		return false;
	}
	//return the computed backslash
	*backslash = focpos2 - focpos1;
	return true;
}