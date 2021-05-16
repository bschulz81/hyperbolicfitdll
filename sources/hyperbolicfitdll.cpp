/// Copyright(c) < 2021 > 
// <Benjamin Schulz> 
// Responsible for:
// The implementation of the library in C++, 
// the development of a parallelized RANSAC inspired algorithm that can removes outlier data in the function "focusposition_Regression",
// the implementation of a repeated median regression in the function "focusposition_Regression" instead of a simple regression,
// the implementation of various known outlier detection methods within the RANSAC 
// (MAD, S, Q, and T, biweight-midvariance estimators), Grubb's test for outliers, Peirce's criterion

// <Stephen King> 
// Responsible for:
// The Algorithm idea for the hyperbolic fit in the internal function "Regression" with the linear regression and least squares comparison, first suggested in
// https://aptforum.com/phpbb/viewtopic.php?p=25998#p25998 

// <Jim Hunt>
// Responsible for:
// Some testing of the library,
// the algorithm idea for "findbackslash_Regression", first suggested in
// https://aptforum.com/phpbb/viewtopic.php?p=26265#p26265 

// The user <cytan>  
// Responsible for:
// The suggestion in https://aptforum.com/phpbb/viewtopic.php?p=26471#p26471 to throw out outlier data by comparison of the error with the Standard-deviation.
// (Note that this old idea is now supplemented by more advanced methods, since the average and standard deviation are not robust.)

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

#include <vector>
#include <chrono>
#include <algorithm>
#include "hyperbolicfitdll.h"
#include <random>  
#include <valarray>
#include <atomic>
#include<unordered_set>
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
inline bool findmodel(valarray<long>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale,valarray<bool>*usedpoint, size_t usedpoints, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method);
inline bool Search_min_error(valarray<long>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression, vector<double>* errs, valarray<bool>* indicestouse);
inline double median(valarray<double>* arr, size_t n, size_t nhalf);
inline double lowmedian(valarray<double>* arr, size_t n);
inline double factorial(size_t n);
inline double peirce(size_t pointnumber, size_t numberofoutliers, size_t fittingparameters);
inline size_t binominal(size_t n, size_t k);
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


inline bool regression(double* ferr, valarray<long>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression, valarray<bool>* usepoint)
{
	*ferr = DBL_MAX;

	valarray<long> x2 = (*x)[*usepoint];
	valarray<double> line_y2 = (*line_y)[*usepoint];
	size_t usedpoints = x2.size();

	valarray<double>line_x2(usedpoints);


	if (use_median_regression)
	{

		valarray<double> stacks2(usedpoints);
		valarray<double> stacks1(usedpoints - 1);

		size_t halfsize = usedpoints / 2;

		valarray<double> stacki1(usedpoints);


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
		double sumy = 0;
		for (size_t n = 0; n < usedpoints; n++)
		{
			sumy += line_y2[n];
		}

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
	*ferr = *ferr /(double) usedpoints;
	return true;
}


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


_inline double factorial(size_t n)
{
	double ret = 1.00;
	for (size_t i = 2; i <= n; ++i)
		ret *= (double)i;
	return ret;
}

_inline double H(double y, size_t nu)
{
	double sum = 0;

	for (size_t j = 1; j <= size_t((nu + 1) / 2 - 1); j++)
	{
		sum += factorial((size_t)j) * factorial(size_t(j - 1.0)) / (pow(4, -((double)j)) * factorial((size_t)2.0 * j)) * pow((1 + y * y / nu), -((double)j));
	}
	return y / (2 * sqrt((double)nu)) * sum;
}
_inline double G(double y, size_t nu)
{
	double sum = 0;
	for (size_t j = 0; j <= nu / 2 - 1; j++)
	{
		sum += factorial((size_t)2.0 * j) / (pow(2.0, 2.0 * j) * pow(factorial(j), 2)) * pow(1 + y * y / nu, -((double)j));
	}
	return sum;
}
inline double cot(double x)
{
	return  cos(x) / sin(x);
}



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

inline double crit(double alpha, size_t N)
{
	double temp = pow(t(alpha / (2.0 * (double)N), N - 2), 2.0);
	return ((double)N - 1.0) / sqrt((double)N) * sqrt(temp / ((double)N - 2.0 + temp));
}

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


inline size_t binominal(size_t n, size_t k)
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

inline double Q_estimator(valarray<double>* err, size_t s)
{
	valarray<double> t1(binominal(s, 2));

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
	size_t k = binominal(h, 2) - 1;

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
inline bool findmodel(valarray<long>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, valarray<bool>* usedpoint, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method)
{
	long this_focpos=0;
	double thiserr=DBL_MAX, thisslope=0.0, thisintercept=0.0;
	if (!Search_min_error(x, line_y, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &this_focpos, use_median_regression, usedpoint))
	{
		return false;
	}
	vector<maybe_inliner> mp;
	size_t pointnumber = (*x).size();
	mp.reserve(pointnumber);

	valarray<double>err(pointnumber);
	valarray<double>err_sqr(pointnumber);

	size_t pointnumberhalf = pointnumber / 2;

	for (size_t p = 0; p < pointnumber; p++)
	{
		double xh = ((double)(*x)[p] - (double)this_focpos);
		xh *= xh;
		double z = fabs(thisslope * xh +thisintercept) - (*line_y)[p];

		err[p] = z;
		err_sqr[p] = z * z;
		if (!(*usedpoint)[p])
		{
			maybe_inliner o;
			o.point = p;
			o.error = z;
			mp.push_back(o);
		}
	}

	if (mp.size() > 0)
	{
		double m = 0, MAD = 0, average = 0, stdev = 0, S = 0, Q = 0, T = 0, biweightmidvariance = 0;

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


		for (size_t j = 0; j < mp.size(); j++)
		{
			bool isoutlier = false;

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

bool focusposition_Regression(vector<long> x, vector<double> y, long* focpos, double* main_error, double* main_slope, double* main_intercept,
	vector<size_t>* indices_of_used_points,
	vector<double>* usedpoints_line_x, vector<double>* usedpoints_line_y, vector<size_t>* indices_of_removedpoints, vector<double>* removedpoints_line_x, vector<double>* removedpoints_line_y,
	double stop_after_seconds, size_t stop_after_numberofiterations_without_improvement, long backslash, double scale, bool use_median_regression,
	size_t maximum_number_of_outliers, outlier_criterion rejection_method, double tolerance)
{
	if (focpos == NULL)
	{
		return false;
	}
	if (x.size() < 4)
	{
		return false;
	}

	size_t pointnumber = x.size();

	size_t minimummodelsize = pointnumber - maximum_number_of_outliers;

	if (maximum_number_of_outliers == 0)
	{
		rejection_method = no_rejection;
	}

	if (rejection_method == no_rejection)
	{
		maximum_number_of_outliers = 0;
		minimummodelsize = pointnumber;
	}

	if (minimummodelsize < 4)
	{
		minimummodelsize = 4;
		maximum_number_of_outliers = pointnumber - minimummodelsize;
	}


	if (rejection_method == tolerance_is_significance_in_Grubbs_test)
	{
		if (tolerance >= 0.9)
		{
			rejection_method = no_rejection;
		}
	}

	long minfocus = LONG_MAX, maxfocus = LONG_MIN;


	valarray<bool> indices(pointnumber);
	valarray<bool> indices2(pointnumber);
	valarray<double>line_yv(pointnumber);

	for (size_t i = 0; i < pointnumber; i++)
	{

		if (i < maximum_number_of_outliers)
		{
			indices[i] = (false);
		}
		else
		{
			indices[i] = (true);
		}
		line_yv[i] = (y[i] * y[i]);
		if (x[i] < minfocus)
		{
			minfocus = x[i];
		}
		if (x[i] > maxfocus)
		{
			maxfocus = x[i];
		}
	}

	if (maxfocus - minfocus == 0)
		return false;


	double  error = DBL_MAX, intercept = 0, slope = 0, additionaldata = 0;


	if (rejection_method == use_peirce_criterion)
	{
		additionaldata = peirce(pointnumber, maximum_number_of_outliers, 3);
	}
	if (rejection_method == tolerance_is_significance_in_Grubbs_test)
	{
		additionaldata = crit(tolerance, pointnumber);
	}


	valarray <long> xv(x.data(), x.size());


	size_t numbercomp = binominal(pointnumber, maximum_number_of_outliers);

#if __cplusplus == 201703L
	std::mutex mtx;
#endif
	


	const size_t number_of_attempts = 705432;



	if (stop_after_numberofiterations_without_improvement < number_of_attempts)
	{
		stop_after_numberofiterations_without_improvement = number_of_attempts;
	}

	if (numbercomp <= number_of_attempts)
	{
		std::unordered_set<std::vector<bool>> helper3;
		vector<valarray<bool>> arr(numbercomp);

		for (size_t i = 0; i < numbercomp; i++)
		{
			arr[i] = indices;
			std::next_permutation(std::begin(indices), std::end(indices));
		}

#if __cplusplus == 201703L
		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method);
			});
		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
			std::vector<bool> helper(pointnumber);
			helper.assign(std::begin(arri), std::end(arri));
			mtx.lock();
			bool b2 = helper3.insert(helper).second;
			mtx.unlock();
			if(b2)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr,&thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
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
			}
		});
#else

		#pragma omp parallel for
		for (long i = 0; i < (long)numbercomp; i++)
		{
			bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method);
		}

		#pragma omp parallel for
		for (long i = 0; i < (long)numbercomp; i++)
		{
			std::vector<bool> helper(pointnumber);
			helper.assign(std::begin(arr[i]), std::end(arr[i]));
			bool b2;
			#pragma omp critical
			{
				b2 = helper3.insert(helper).second;
			}

			if(b2)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept,&thisslope, &thisfocpos, use_median_regression, &arr[i]);
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
							indices2 = arr[i];
						}
					}
				}
			}
		}
#endif
	}
	else if (numbercomp <= number_of_attempts * 100)
	{
		size_t p =(size_t) numbercomp/number_of_attempts;
		
		for (size_t o = 0; o < p; o++)
		{
			std::unordered_set<std::vector<bool>> helper3;
			vector<valarray<bool>> arr(number_of_attempts);

			for (size_t i = 0; i < number_of_attempts; i++)
			{
				arr[i] = indices;
				std::next_permutation(std::begin(indices), std::end(indices));
			}

#if __cplusplus == 201703L
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method);
				});
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					std::vector<bool> helper(pointnumber);
					helper.assign(std::begin(arri), std::end(arri));
					mtx.lock();
					bool b2 = helper3.insert(helper).second;
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
								indices2 = arri;
							}
							mtx.unlock();
						}
					}
				});
#else	
			#pragma omp parallel for
			for (long i = 0; i < (long)number_of_attempts; i++)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method);
			}

			#pragma omp parallel for
			for (long i = 0; i < (long)number_of_attempts; i++)
			{
				std::vector<bool> helper(pointnumber);
				helper.assign(std::begin(arr[i]), std::end(arr[i]));
				bool b2;
				#pragma omp critical
				{
					b2 = helper3.insert(helper).second;
				}

				if (b2)
				{
					double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
					long thisfocpos = 0;
					bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr[i]);
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
								indices2 = arr[i];
							}
						}
					}
				}
			}
#endif
		}

		std::unordered_set<std::vector<bool>> helper3;
		vector<valarray<bool>> arr((numbercomp % number_of_attempts)+1);

		for (size_t i = 0; i <= (numbercomp % number_of_attempts); i++)
		{
			arr[i] = indices;
			std::next_permutation(std::begin(indices), std::end(indices));
		}
#if __cplusplus == 201703L
		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method);
			});
		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				std::vector<bool> helper(pointnumber);
				helper.assign(std::begin(arri), std::end(arri));
				mtx.lock();
				bool b2 = helper3.insert(helper).second;
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
							indices2 = arri;
						}
						mtx.unlock();
					}
				}
			});
#else
	
		#pragma omp parallel for
		for (long i = 0; i <= (long)(numbercomp % number_of_attempts); i++)
		{
			bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method);
		}

		#pragma omp parallel for
		for (long i = 0; i <= (long)(numbercomp % number_of_attempts); i++)
		{
			std::vector<bool> helper(pointnumber);
			helper.assign(std::begin(arr[i]), std::end(arr[i]));
			bool b2;
			#pragma omp critical
			{
				b2 = helper3.insert(helper).second;
			}

			if (b2)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr[i]);
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
							indices2 = arr[i];
						}
					}
				}
			}
		}
#endif
	}

	
	else
	{
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
			std::unordered_set<std::vector<bool>> helper3;
			valarray<valarray<bool>> arr(number_of_attempts);
			for (size_t i = 0; i < number_of_attempts; i++)
			{
				shuffle(std::begin(indices), std::end(indices), urng);
				arr[i] = indices;
			}

#if __cplusplus == 201703L
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method);
				});
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					counter1++;
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
				bool b1 = findmodel(&xv, &line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method);
			}

			#pragma omp parallel for
			for (long i = 0; i < number_of_attempts; i++)
			{
				#pragma omp atomic
				counter1++;

				std::vector<bool> helper;
				helper.assign(std::begin(arr[i]), std::end(arr[i]));
				bool b2;
				#pragma omp critical
				{
					b2 = helper3.insert(helper).second;
				}

				if (b2)
				{
					double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
					long thisfocpos = 0;
					bool b3 = Search_min_error(&xv, &line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr[i]);
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
								indices2 = arr[i];
								counter1 = 0;
							}
						}
					}
				}
			}
#endif

			if (counter1 >= stop_after_numberofiterations_without_improvement)
			{
				break;
			}
			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			seconds = elapsed_seconds.count();
		} while (seconds < fabs(stop_after_seconds));
	}
	

	if (error == DBL_MAX)
	{
		return false;
	}

	*focpos += backslash;

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

			if(removedpoints_line_y != NULL)
			{
				double w = y[k];
				(*removedpoints_line_y).push_back(w * w);
			}
		}

	}


	return true;
}


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
	if (focusposition_Regression(x1, y1, &focpos1, main_error1, main_slope1, main_intercept1, indicesofusedpoints1, used_points1_line_x, used_points1_line_y, indicesofremovedpoints1, removedpoints1_line_x, removedpoints1_line_y, stop_after_seconds,
		stop_after_numberofiterations_without_improvement, 0, scale, use_median_regression, maximum_number_of_outliers, rejection_method, tolerance) == false)
	{
		return false;
	}

	if (focusposition_Regression(x2, y2, &focpos2, main_error2, main_slope2, main_intercept2, indicesofusedpoints2, used_points2_line_x, used_points2_line_y, indicesofremovedpoints2, removedpoints2_line_x, removedpoints2_line_y, stop_after_seconds,
		stop_after_numberofiterations_without_improvement, 0, scale, use_median_regression, maximum_number_of_outliers, rejection_method, tolerance) == false)
		*backslash = focpos2 - focpos1;
	return true;
}



