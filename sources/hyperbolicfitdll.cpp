// Copyright(c) < 2021 > 

// <Benjamin Schulz> 
// Responsible for:
// The implementation of the library in C++, see https://aptforum.com/phpbb/viewtopic.php?p=26587#p26587),
// The development of a RANSAC inspired algorithm to remove outliers in the function "focusposition_Regression", 
// The implementation of a repeated median regression in the function "focusposition_Regression" instead of a simple regression,
// The implementation of various known outlier detection detection methods within the RANSAC (MAD, S, Q, and T estimators), Grubb's test for outliers, Peirce's criterion.

// <Stephen King> 
// Responsible for:
// The Algorithm idea for the hyperbolic fit in the internal function "Regression" with the linear regression and least squares comparison, first suggested in  https://aptforum.com/phpbb/viewtopic.php?p=25998#p25998).

// <Jim Hunt>
// Responsible for:
// Some testing of the library,
// the algorithm idea for "findbackslash_Regression", first suggested in https://aptforum.com/phpbb/viewtopic.php?p=26265#p26265).

// The user <cytan> in https://aptforum.com/phpbb/viewtopic.php?p=26471#p26471 
// Responsible for:
// The suggestion to throw out outlier data in "focusposition_Regression" by comparison of the error with the Standard-deviation.
// Note that this idea is now supplemented by more advanced methods, since the average and standard deviation are not robust.

// The library makes use of an algorithm for student's distribution, which can be found at
// Smiley W. Cheng, James C. Fu, Statistics & Probability Letters 1 (1983), 223-227

// The library also makes use of Peirce's outlier test. An algorithm for this method was developed in
// Gould, B. A, Astronomical Journal, vol. 4, iss. 83, p. 81 - 87 (1855).

// The library also has the possibility to use MAD, S and Q estimators. 
// These estimators are extensively described in Peter J. Rousseeuw, Christophe Croux, Alternatives to the Median-Absolute Deviation
// J. of the Amer. Statistical Assoc. (Theory and Methods), 88 (1993),p. 1273,
// Christophe Croux and Peter J.Rousseeuw, Time-effcient algorithms for two highly robust estimators of scale, 
// In: Dodge Y., Whittaker J. (eds) Computational Statistics. Physica, Heidelberg, https ://doi.org/10.1007/978-3-662-26811-7_58


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
struct maybe_inliner
{
	size_t point;
	double error;
};

inline void stdeviation(vector<double>* errs, double* stdev, double* average);
inline bool Ransac_regression(vector<long>* x, vector<double>* y, vector<double>* line_y, size_t pointnumber, long minfocus, long maxfocus, double scale,
	vector<bool>usedpoint,size_t usedpoints, double tolerance, long* this_focpos, double* thiserr,	double* thisslope, double* thisintercept, vector<size_t>* usedindices, 
	vector<size_t>* removedindices, double additionaldata, bool use_median_regression, outlier_criterion rejection_method);
inline bool regression(double* ferr, vector<long>* x, vector<double>* y, vector<double>* line_y, size_t pointnumber, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression,size_t usedpoints);
inline bool Search_min_error(vector<long>* x, vector<double>* y, vector<double>* line_y, size_t pointnumber, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression, vector<double>* errs, vector<bool>* indicestouse,size_t usedpoints);
inline bool interceptslope(vector<double>* x, vector<double>* y, double* fintercept, double* fslope, bool use_median_regression, vector<bool>* usepoint, size_t usedpoints);
inline double median(vector<double>* arr);
inline double median(vector<double>* arr, size_t n, size_t nhalf);
inline double lowmedian(vector<double>*arr);
inline double factorial(size_t n);
inline double peirce(size_t pointnumber, size_t numberofoutliers, size_t fittingparameters);
inline size_t binominal(size_t n, size_t k);
inline double crit(double alpha, size_t N);
inline double t(double alpha, size_t nu);
inline double cot(double x);
inline double G(double y, size_t nu);
inline double H(double y, size_t nu);
inline double Q_estimator(vector<double>* err);
inline double S_estimator(vector<double>* err);
inline double MAD_estimator(vector<double>* err, double* m);
inline double T_estimator(vector<double>* err);


inline void stdeviation(vector<double>* errs, double* stdev, double* average)
{
	double sum = 0, devi = 0, t = 0;
	size_t s = (long)(*errs).size();
	for (size_t i = 0; i <s; i++)
	{
		sum += (*errs)[i];
	}
	*average = sum /(double) s;

	if (stdev != NULL)
	{
		for (size_t i = 0; i <s; i++)
		{
			t = ((*errs)[i] - *average);
			t *= t;
			devi += t;
		}
		*stdev = sqrt(devi / s);
	}
}




inline double median(vector<double>* arr)
{
	size_t n = (*arr).size();
	size_t nhalf = n / 2;
	nth_element((*arr).begin(), (*arr).begin() + nhalf, (*arr).end());
	double  med = (*arr)[nhalf];
	if (n % 2 == 0) 
	{ 
		auto max_it = max_element((*arr).begin(), (*arr).begin() + nhalf);
		med = (*max_it + med) / 2.0;
	}
	return med;
}
inline double median(vector<double>* arr,size_t n, size_t nhalf)
{
	nth_element((*arr).begin(), (*arr).begin() + nhalf, (*arr).end());
	double  med = (*arr)[nhalf];
	if (n % 2 == 0)
	{
		auto max_it = max_element((*arr).begin(), (*arr).begin() + nhalf);
		med = (*max_it + med) / 2.0;
	}
	return med;
}


inline double lowmedian(vector<double>* arr)
{
	size_t n = (*arr).size();
	size_t m = (size_t)(floor(((double)n + 1.0) / 2.0) - 1.0);
	nth_element((*arr).begin(), (*arr).begin() + m, (*arr).end());
	return (double)(*arr)[m];
}



inline bool interceptslope(vector<double>* x, vector<double>* y, double* fintercept, double* fslope, bool use_median_regression,vector<bool>*usepoint,size_t usedpoints)
{
	size_t s = (*x).size();
	if (use_median_regression)
	{
		vector<double> stacks2;
		stacks2.resize(usedpoints);
		vector<double> stacks1;
		size_t size2 = usedpoints - 1;
		stacks1.resize(size2);

		size_t k = 0,halfsize=size2/2;
		for (size_t i = 0; i < s; i++)
		{
			size_t q = 0;
			if ((*usepoint)[i])
			{
				for (size_t j = 0; j < s; j++)
				{
					if (i != j)
					{
						if ((*usepoint)[j])
						{
							stacks1[q] = ((*y)[j] - (*y)[i]) / ((*x)[j] - (*x)[i]);
							q++;
						}
					}
				}
				stacks2[k] =median(&stacks1,size2,halfsize);
				k++;
			}
		}
		halfsize = usedpoints / 2;
		*fslope = median(&stacks2,usedpoints,halfsize);

		vector<double> stacki1;
		stacki1.resize(usedpoints);

		k = 0;

		for (size_t i = 0; i < s; i++)
		{
			if ((*usepoint)[i])
			{
				stacki1[k] = ((*y)[i] - *fslope * (*x)[i]);
				k++;
			}
		}
		*fintercept = median(&stacki1,usedpoints,halfsize);
	}
	else
	{
		double sumx = 0, sumy = 0, sumxy = 0, sumxx = 0;
		for (size_t i = 0; i < s; i++)
		{
			if ((*usepoint)[i])
			{
				sumx += (double)(*x)[i];
				sumy += (double)(*y)[i];
				sumxy += ((double)(*x)[i]) * ((double)(*y)[i]);
				sumxx += ((double)(*x)[i]) * ((double)(*x)[i]);
			}
		}

		double xaverage = sumx / (double)usedpoints;
		double yaverage = sumy / (double)usedpoints;

		double t = sumxx - sumx * xaverage;
		if (t == 0)
		{
			return false;
		}
		*fslope = (sumxy - sumx * yaverage) / t;
		*fintercept = yaverage - (*fslope) * xaverage;
	}
	return true;
}

inline bool regression(double* ferr, vector<long>* x, vector<double>* y, vector<double>* line_y, size_t pointnumber, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression,vector<bool>*usepoint,size_t usedpoints)
{
	
	*ferr = DBL_MAX;
	vector<double>line_x;
	line_x.resize(pointnumber);


	for (long h = minfocus; h <= maxfocus; h++)
	{
		for (size_t n = 0; n < pointnumber; n++)
		{
			if ((*usepoint)[n])
			{
				double k = (double)(*x)[n] - (double)h;
				line_x[n] = k * k;
			}
		}
		double thiserr=0, thisslope=0, thisintercept=0;
		if (interceptslope(&line_x, line_y, &thisintercept, &thisslope, use_median_regression,usepoint,usedpoints) == false)
		{
			return false;
		}

		for (size_t n = 0; n < pointnumber; n++)
		{
			if ((*usepoint)[n])
			{
				double k = sqrt(fabs(thisslope * (double)line_x[n] + thisintercept)) - (*y)[n];
				k *= k;
				thiserr += k;
			}
		}
		if (thiserr < *ferr)
		{
			*fintercept = thisintercept;
			*fslope = thisslope;
			*ferr = thiserr;
			*focpos = h;
		}
	}
	
	return true;
}


inline bool Search_min_error(vector<long>* x, vector<double>* y,vector<double>*line_y, size_t pointnumber, long minfocus, long maxfocus,double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression,vector<bool>*indicestouse,size_t usedpoints)
{

	regression(err, x, y, line_y, pointnumber, minfocus, maxfocus, focpos, fintercept, fslope, use_median_regression,indicestouse,usedpoints);

	if ((*focpos == maxfocus)&& (scale >1))
	{
		
		const long middle=lround(((double)maxfocus +(double) minfocus) / 2.0);
		minfocus = maxfocus + 1;
		maxfocus = (long)((double)middle + ((double)maxfocus -(double) middle) * fabs(scale));
		long thisfocpos;
		double thiserr, thisslope, thisintercept;
		regression(&thiserr,x, y, line_y, pointnumber, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope, use_median_regression, indicestouse,usedpoints);
	
		if ((thiserr < *err))
		{
			*err = thiserr;
			*fintercept = thisintercept;
			*fslope = thisslope;
			*focpos = thisfocpos;
		}
	}
	else if ((*focpos == minfocus )&& (scale >1) )
	{
		const long middle = lround(((double)maxfocus + (double)minfocus) / 2.0);
		maxfocus = minfocus - 1;
		minfocus = (long)((double)middle - ((double)middle-(double)minfocus) * fabs(scale));
		long thisfocpos;
		double thiserr, thisslope, thisintercept;
		regression(&thiserr,x, y, line_y, pointnumber, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope,use_median_regression, indicestouse,usedpoints);
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


inline double factorial(size_t n)
{
	double ret = 1.00;
	for (size_t i = 2; i <= n; ++i)
		ret *= (double)i;
	return ret;
}

inline double H(double y, size_t nu)
{
	double sum = 0;

	for (size_t j = 1; j <=size_t((nu + 1) / 2 - 1); j++)
	{
		sum += factorial((size_t)j) * factorial(size_t(j - 1.0)) / (pow(4, -((double)j)) * factorial((size_t)2.0 * j)) * pow((1 + y * y / nu), -((double)j));
	}
	return y / (2 * sqrt((double)nu)) * sum;
}
inline double G(double y, size_t nu)
{
	double sum = 0;
	for (size_t j = 0; j <=  nu / 2 - 1; j++)
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
			xa = 1.0 + (diff1 - p) / n * (1.0 -lambda * lambda);
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
		size_t prod1 = 1, prod2 = n;
		for (size_t i = 2; i <= k; i++)
		{
			prod1 *= i;
			prod2 *= (n + 1 - i);
		}
		return(size_t)(prod2 / prod1);
	}
}

inline double Q_estimator(vector<double>* err)
{
	size_t s = (*err).size();
	vector<double> t1;
	t1.resize(binominal(s, 2));
	size_t i = 0;
	for (size_t w = 0; w < s - 1; w++)
	{
		for (size_t k = w + 1; k < s; k++)
		{
			t1[i]=(fabs((*err)[w] - (*err)[k]));
			i++;
		}
	}
	size_t h = s / 2 + 1;
	size_t k = binominal(h, 2) - 1;

	std::nth_element(t1.begin(), t1.begin() + k, t1.end());

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
			cc =s / (s+ 3.8);
		}
	}
	return cc * 2.2219 * t1[k];



}

inline double S_estimator(vector<double>* err)
{
	vector<double>m1;
	
	
	size_t s = (*err).size();
	m1.resize(s);

	vector<double>t1;
	t1.resize(s - 1);

	for (size_t i = 0; i < s; i++)
	{
		size_t q = 0;
		for (size_t k = 0; k < s; k++)
		{
			if (i != k)
			{
				t1[q]=(fabs((*err)[i] - (*err)[k]));
				q++;
			}
		}
		m1[i]=(lowmedian(&t1));
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

	return (c * 1.1926 * lowmedian(&m1));
}

inline double MAD_estimator(vector<double>* err,double *m)
{
	*m = median(err);
	vector<double>m1;

	size_t s = (*err).size();
	m1.resize(s);
	for (size_t i = 0; i < s; i++)
	{
		m1[i]=(fabs((*err)[i] - *m));
	}
	double cn[] = { 0, 0, 1.196 ,1.495 ,1.363, 1.206, 1.200, 1.140, 1.129,1.107 };
	double c = 1.0;
	if (s <= 9)
	{
		c = cn[s];
	}
	else
	{
	    c =s/ (s - 0.8);
	}

	return  c * 1.4826 * median(&m1);
}

inline double T_estimator(vector<double>* err)
{
	
	vector<double>med1;
	size_t s = (*err).size();
	med1.resize(s);
	vector<double>arr;
	arr.resize(s - 1);

	for (size_t i = 0; i <s; i++)
	{
		size_t q = 0;
	
		for (size_t j = 0; j <s; j++)
		{
			if (i != j)
			{
				arr[q]=(fabs((*err)[i] - (*err)[j]));
				q++;
			}
		}
		med1[i]=(median(&arr));
	}
	size_t h = s / 2 + 1;
	double w = 0;
	sort(med1.begin(), med1.end());
	for (size_t i = 1; i <= h; i++)
	{
		w += med1[i - 1];
	}
	return  (1.38 /((double) h))*w;
}
inline bool Ransac_regression(vector<long>* x, vector<double>* y, vector<double>* line_y, size_t pointnumber, long minfocus, long maxfocus, double scale, vector<bool>usedpoint,size_t usedpoints, double tolerance, long* this_focpos, double* thiserr,
	double* thisslope, double* thisintercept, vector<size_t>* usedindices, vector<size_t>* removedindices,double additionaldata, bool use_median_regression, outlier_criterion rejection_method)
{

	if (Search_min_error(x, y, line_y, pointnumber, minfocus, maxfocus, scale, thiserr, thisintercept, thisslope, this_focpos, use_median_regression, &usedpoint,usedpoints) == false)
	{
		return false;
	}
	vector<maybe_inliner> mp;
	vector<double>err;
	vector<double>err_sqr;


	err.resize(pointnumber);
	err_sqr.resize(pointnumber);

	mp.reserve(pointnumber);

	if (removedindices != NULL)
	{
		(*removedindices).clear();
		(*removedindices).reserve(pointnumber - 4);
	}
	if (usedindices != NULL)
	{
		(*usedindices).clear();
		(*usedindices).reserve(pointnumber);
	}


	for (size_t p = 0; p < (*x).size(); p++)
	{
		double xh = ((double)(*x)[p] - (double)*this_focpos);
		xh *= xh;
		double z =fabs(*thisslope * xh + *thisintercept) - (*y)[p]*(*y)[p];

		err[p]=z;
		err_sqr[p]=z * z;
		if (usedpoint[p])
		{
			if (usedindices != NULL)
			{
				(*usedindices).push_back(p);
			}
		}
		else
		{
			maybe_inliner o;
			o.point = p;
			o.error = z;
			mp.push_back(o);
		}
	}

	if (mp.size() > 0)
	{
		double m = 0, MAD = 0, average = 0, stdev = 0, S = 0, Q = 0,T=0;

		switch (rejection_method)
		{
			case tolerance_multiplies_standard_deviation_of_error:
			{
				stdeviation(&err, &stdev, &average);
				break;
			}
			case tolerance_is_significance_in_Grubbs_test:
			{
				stdeviation(&err, &stdev, &average);
				break;
			}
			case tolerance_is_decision_in_MAD_ESTIMATION:
			{
				MAD = MAD_estimator(&err, &m);
				break;
			}
			case tolerance_is_decision_in_S_ESTIMATION:
			{
				S = S_estimator(&err);
				m = median(&err);
				break;
			}
			case tolerance_is_decision_in_Q_ESTIMATION:
			{
				m = median(&err);
				Q = Q_estimator(&err);
				break;
			}
			case tolerance_is_decision_in_T_ESTIMATION:
			{
				m = median(&err);
				T = T_estimator(&err);
				break;
			}
			case use_peirce_criterion:
			{
				stdeviation(&err_sqr, NULL, &average);
				break;
			}
		}

		for (size_t j = 0; j < mp.size(); j++)
		{
			bool isoutlier=false;
		
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

					double G = fabs(average - mp[j].error) / stdev;
					if (G > additionaldata)
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
				case tolerance_is_decision_in_S_ESTIMATION:
				{
					double G = (fabs(mp[j].error - m) / S);
					if (G> tolerance)
					{
						isoutlier = true;
					}
					break;
				}

				case tolerance_is_decision_in_Q_ESTIMATION:
				{
					double G= fabs((mp[j].error - m) / Q);
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

			if(isoutlier)
			{
				usedpoint[mp[j].point] = false;
				if (removedindices != NULL)
				{
					(*removedindices).push_back(mp[j].point);
				}
			}
			else
			{
				usedpoints++;
				usedpoint[mp[j].point] = true;
				if (usedindices != NULL)
				{
					(*usedindices).push_back(mp[j].point);
				}
			}
		}
	}
	if (removedindices != NULL)
	{
		(*removedindices).shrink_to_fit();
	}
	if (usedindices != NULL)
	{
		(*usedindices).shrink_to_fit();
	}

	if (Search_min_error(x, y, line_y, pointnumber, minfocus, maxfocus, scale, thiserr, thisintercept, thisslope, this_focpos, use_median_regression, &usedpoint,usedpoints) == false)
	{
		return false;
	}
	return true;
}

bool focusposition_Regression(vector<long> x, vector<double> y, long* focpos, double* main_error, double* main_slope, double* main_intercept,
	vector<size_t>* indices_of_used_points,
	vector<double>* usedpoints_line_x, vector<double>* usedpoints_line_y, vector<size_t>* indices_of_removedpoints, vector<double>* removedpoints_line_x, vector<double>* removedpoints_line_y,
	double stop_after_seconds , size_t stop_after_numberofiterations_without_improvement, long backslash, double scale, bool use_median_regression,
	size_t maximum_number_of_outliers, outlier_criterion rejection_method, double tolerance)
{
	if(focpos==NULL )
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


	vector<bool> indices;
	
	vector<double>line_y;

	line_y.resize(pointnumber);
	indices.resize(pointnumber);
	for (size_t i = 0; i < pointnumber; i++)
	{

		if (i < maximum_number_of_outliers)
		{
			indices[i]=(false);
		}
		else 
		{
			indices[i]=(true);
		}
		line_y[i]=(y[i] * y[i]);
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

	vector<size_t> removedindices, usedindices, thisremovedindices, thisusedindices;

	double thiserr = DBL_MAX, thisslope=0, thisintercept=0, error=DBL_MAX,intercept=0,slope=0, additionaldata = 0;
	long thisfocpos = 0;

	if (rejection_method == use_peirce_criterion)
	{
		additionaldata = peirce(pointnumber, maximum_number_of_outliers, 3);
	}
	if (rejection_method==tolerance_is_significance_in_Grubbs_test)
	{
		additionaldata = crit(tolerance, pointnumber);
	}

	double k;
	if (binominal(pointnumber,maximum_number_of_outliers)<=(size_t) 184756)
	{
		do
		{
			Ransac_regression(&x, &y, &line_y, pointnumber, minfocus, maxfocus, scale, indices, minimummodelsize, tolerance, &thisfocpos, &thiserr, &thisslope, &thisintercept, &thisusedindices, &thisremovedindices, additionaldata, use_median_regression, rejection_method);
			k = thiserr / thisusedindices.size();
			if ( k< error)
			{
				error = k;
				*focpos = thisfocpos;
				slope = thisslope;
				intercept = thisintercept;
				usedindices = thisusedindices;
				removedindices = thisremovedindices;
			}
		} while (std::next_permutation(indices.begin(), indices.end()));
	}
	else
	{
		double seconds = 0;
		size_t counter1 = 0;
		std::random_device rng;
		std::mt19937 urng(rng());
		auto start = std::chrono::steady_clock::now();
		do
		{
			Ransac_regression(&x, &y,&line_y, pointnumber,minfocus,maxfocus,scale, indices, minimummodelsize, tolerance, &thisfocpos, &thiserr, &thisslope, &thisintercept,&thisusedindices, &thisremovedindices, additionaldata, use_median_regression, rejection_method);
			k = thiserr / thisusedindices.size();
			if (k < error)
			{
				error = k;
				*focpos = thisfocpos;
				slope = thisslope;
				intercept = thisintercept;
				usedindices = thisusedindices;
				removedindices = thisremovedindices;
				counter1 = 0;
			}
			shuffle(indices.begin(), indices.end(), urng);
			if (counter1 == stop_after_numberofiterations_without_improvement)
			{
				break;
			}
			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			seconds = elapsed_seconds.count();		
			counter1++;
		} while (seconds < fabs(stop_after_seconds));
	}


	*focpos += backslash;

	if (main_error != NULL)
	{
		*main_error =error;
	}

	if (main_slope != NULL)
	{
		*main_slope =slope;
	}

	if (main_intercept != NULL)
	{
		*main_intercept =intercept;
	}

	if (indices_of_removedpoints != NULL)
	{
		*indices_of_removedpoints = removedindices;
	}

	if (indices_of_used_points != NULL)
	{
		*indices_of_used_points = usedindices;
	}

	if((removedpoints_line_x != NULL) && (removedpoints_line_y != NULL))
	{
		(*removedpoints_line_x).clear();
		(*removedpoints_line_y).clear();
		(*removedpoints_line_x).resize(removedindices.size());
		(*removedpoints_line_y).resize(removedindices.size());
		for (size_t c = 0; c < removedindices.size(); c++)
		{
			k =(double) x[removedindices[c]] -(double) *focpos;
			k *= k;
			(*removedpoints_line_x)[c]=(k);
			double w = y[removedindices[c]];
			(*removedpoints_line_y)[c]=(w * w);
		}
	}
	if ((usedpoints_line_x != NULL) && (usedpoints_line_y != NULL))
	{
		(*usedpoints_line_x).clear();
		(*usedpoints_line_y).clear();
		(*usedpoints_line_x).resize(usedindices.size());
		(*usedpoints_line_y).resize(usedindices.size());
		for (size_t c = 0; c < usedindices.size(); c++)
		{
			k = (double)x[usedindices[c]] -(double) *focpos;
			k *= k;
			(*usedpoints_line_x)[c]=(k);
			double w = y[usedindices[c]];
			(*usedpoints_line_y)[c]=(w*w);
		}
	}


	return true;
}


bool findbackslash_Regression(long* backslash,
	vector<long> x1, vector<double> y1, vector<long> x2, vector<double> y2, double* main_error1 , double* main_slope1, double* main_intercept1, vector<size_t>*indicesofusedpoints1,
	vector<double>*used_points1_line_x , vector<double>*used_points1_line_y , vector<size_t>*indicesofremovedpoints1, vector<double>*removedpoints1_line_x, vector<double>*removedpoints1_line_y ,
	double* main_error2 , double* main_slope2 , double* main_intercept2 , vector<size_t>*indicesofusedpoints2 ,
	vector<double>*used_points2_line_x , vector<double>*used_points2_line_y , vector<size_t>*indicesofremovedpoints2, vector<double>*removedpoints2_line_x, vector<double>*removedpoints2_line_y ,
	double stop_after_seconds, size_t stop_after_numberofiterations_without_improvement, double scale, bool use_median_regression,
	size_t maximum_number_of_outliers, outlier_criterion rejection_method, double tolerance)
{
	if (backslash == NULL)
	{
		return false;
	}

	long focpos1 = 0, focpos2 = 0;
	if (focusposition_Regression(x1, y1, &focpos1, main_error1, main_slope1, main_intercept1,indicesofusedpoints1, used_points1_line_x, used_points1_line_y, indicesofremovedpoints1, removedpoints1_line_x, removedpoints1_line_y, stop_after_seconds,
		stop_after_numberofiterations_without_improvement, 0, scale,use_median_regression,maximum_number_of_outliers,rejection_method,tolerance) == false)
	{
		return false;
	}

	if (focusposition_Regression(x2, y2, &focpos2, main_error2, main_slope2, main_intercept2,indicesofusedpoints2, used_points2_line_x, used_points2_line_y, indicesofremovedpoints2, removedpoints2_line_x, removedpoints2_line_y, stop_after_seconds,
		stop_after_numberofiterations_without_improvement,0, scale, use_median_regression, maximum_number_of_outliers, rejection_method, tolerance) == false)
		*backslash = focpos2 - focpos1;
	return true;
}


