// This file implements the algorithms and procedures described inMonthly Notices of the Royal Astronomical Society,
// Volume 511, Issue 2, April 2022, Pages 2008–2020, https://doi.org/10.1093/mnras/stac189 
// ( for a preprint, see https://arxiv.org/abs/2201.12466 ).

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



//


// The library makes use of an imaging analysis algorithm that was developed by C. Y. Tan (first author) and B. Schulz (second author).

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
// of this software and associated documentation files(the "Software"), to deal
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

#include "focusinterpolation.h"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <random>  
#include <valarray>
#include <atomic>
#include <unordered_set>
#include <complex>
#include <cmath>
#include <omp.h>
using namespace cv;
#if __cplusplus == 201703L
	#include <mutex>
	#include <execution>
#endif


#ifdef GNUCOMPILER
#define _inline inline
#endif

#ifdef CLANGCOMPILER
#define _inline inline
#endif

struct maybe_inliner
{
	size_t point;
	double error;
};

class Vector
{
private:
	std::valarray<double> m;
public:
	inline Vector(size_t);
	inline Vector operator+(Vector);
	inline valarray<double> array();
	inline Vector operator-(Vector);
	inline double operator*(Vector);
	inline Vector operator*(double);
	inline Vector operator/(double);
	inline double& operator()(size_t);
	inline size_t size() const;
	inline void push_back(double t);
	inline void resize(size_t t);
	inline void printvector();
};

class Matrix
{
private:
	size_t r;
	size_t c;
	std::valarray<double> m;
public:
	inline Matrix(size_t, size_t);
	inline Matrix operator+(Matrix);
	inline Matrix operator-(Matrix);
	inline Matrix operator*(Matrix);
	inline Vector operator*(Vector);
	inline Matrix Transpose();
	inline Matrix operator*(double);
	inline Matrix operator/(double);
	inline double& operator()(size_t, size_t);
	inline size_t Rows() const;
	inline size_t Columns() const;
	inline Vector Gaussian_algorithm(Matrix* m, Vector* v);
	inline Matrix Identity(size_t rows, size_t columns);
	inline Matrix Diagonal();
	inline void printmatrix();
	inline void resize(size_t rows,size_t columns);
};


inline Matrix Jacobian(Vector* X, Vector* beta,bool withbe,bool withtheta);
inline Vector f(Vector* X, Vector* beta,bool withbe,bool withtheta);
inline Vector directionalderivative(Vector* X, Vector*fxbeta, Vector* beta, Vector* delta, Matrix* J,bool withbe,bool withtheta);
inline double error(Vector* Y, Vector* X, Vector* beta,bool withbe, bool withtheta);
inline double error(Vector* Y, Vector* X, Vector* beta, bool withbe, bool withtheta, outlier_criterion method, double tolerance, double additionaldata,valarray<bool>*indices);
inline double error(valarray<double>* errs, outlier_criterion method, double tolerance, double additionaldata,valarray<bool>*indicestouse);
inline void computedelta1(Vector* delta, Matrix G, Vector v, double lambda);
inline void computedelta2(Vector* delta, Matrix G, Vector v, double lambda);
inline void stdeviation(valarray<double>* errs, double* stdev, double* average, size_t s);
inline void hyperbola_linear_regression(double* ferr, valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression, valarray<bool>* usepoint);
inline void hyperbola_linear_regression(double* ferr, valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression,outlier_criterion method, double tolerance,double additionaldata, valarray<bool>* usepoint);

inline void hyperbola_regression_extendsearcharea(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression, valarray<bool>* indicestouse);
inline void hyperbola_regression_extendsearcharea(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression,outlier_criterion method, double tolerance, double additionaldata, valarray<bool>* indicestouse);

inline void hyperbola_Least_trimmed_squares(double* error, long* focpos, double* slope, double* intercept, valarray<bool>* returnindices, size_t pointnumber, size_t pointnumberhalf, size_t maximum_number_of_outliers, size_t stop_after_numberofiterations_without_improvement, double stop_after_seconds, valarray<bool>* indices, valarray<double>* xv, valarray<double>* line_yv, long minfocus, long maxfocus, double scale, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method);
inline void hyperbola_findmodel(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, valarray<bool>* usedpoint, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method, size_t pointnumber, size_t pointnumberhalf);
inline double median(valarray<double> *arr, size_t n, size_t nhalf);
inline double median(valarray<float>* arr, size_t n, size_t nhalf);
inline double lowmedian(valarray<double> *arr, size_t n);
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
inline void powerfunction_linear_regression(double* ferr, valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression, valarray<bool>* usepoint);
inline void powerfunction_linear_regression(double* ferr, valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression,outlier_criterion method, double tolerance, double additionaldata, valarray<bool>* usepoint);

inline void powerfunction_findmodel(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, valarray<bool>* usedpoint, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method, size_t pointnumber, size_t pointnumberhalf);
inline void powerfunction_findmodel(valarray<double>* x, valarray<double>* line_y, double alpha, double gamma,long focpos, long minfocus, long maxfocus, double scale, valarray<bool>* usedpoint, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method, size_t pointnumber, size_t pointnumberhalf, size_t stop_after_iterations, double stop_after_seconds,bool withbe, bool withtheta);
inline void powerfunction_regression_extendsearcharea(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression, valarray<bool>* indicestouse);
inline void powerfunction_regression_extendsearcharea(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression,outlier_criterion method,double tolerance, double additionaldata, valarray<bool>* indicestouse);

inline void powerfunction_Least_trimmed_squares(double* error, long* focpos, double* slope, double* intercept, valarray<bool>* returnindices, size_t pointnumber, size_t pointnumberhalf, size_t maximum_number_of_outliers, size_t stop_after_numberofiterations_without_improvement, double stop_after_seconds, valarray<bool>* indices, valarray<double>* xv, valarray<double>* line_yv, long minfocus, long maxfocus, double scale, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method);
inline void powerfunction_nonlinear_regression(double* err, valarray<double>* x, valarray<double>* y, double* alpha, double* gamma, double* b, double* theta, long* focpos, outlier_criterion method, double tolerance, double additionaldata, valarray<bool>* usedpoint, long minfocus, long maxfocus, double scale, size_t stop_after_iterations, double stop_after_seconds, bool withbe, bool withtheta);

inline void powerfunction_nonlinear_regression(double* err, valarray<double>* x, valarray<double>* y, double* alpha, double* gamma, double *b, double* theta, long* focpos, valarray<bool>* usedpoint, long minfocus, long maxfocus,  double scale, size_t stop_after_iterations, double stop_after_seconds,bool withbe,bool withtheta);



inline Vector::Vector(size_t s) {
	m.resize(s,0.0);

}
inline void Vector::resize(size_t t)
{
	m.resize(t);
}

inline valarray<double> Vector::array()
{
	return m;
}
inline Vector Vector::operator+(Vector B) {
	Vector sum(m.size());

#pragma omp parallel for
	for (long i = 0; i < m.size(); i++)
	{
		sum(i) = m[i] + B(i);
	}
	return sum;
}
inline Vector Vector::operator-(Vector B) {
	Vector sum(m.size());

#pragma omp parallel for
	for (long i = 0; i < m.size(); i++)
	{
		sum(i) = m[i] - B(i);
	}
	return sum;
}
inline double Vector::operator*(Vector B) {
	double sum = 0;

#pragma omp parallel for
	for (long i = 0; i < m.size(); i++)
	{
		sum += m[i] * B(i);
	}
	return sum;
}
inline Vector Vector::operator*(double B) {
	Vector sum(m.size());
#pragma omp parallel for
	for (long i = 0; i < m.size(); i++)
	{
		sum(i) = m[i] * B;
	}
	return sum;
}
inline Vector Vector::operator/(double B) {
	Vector sum(this->size());
	if (B != 0)
	{
#pragma omp parallel for
		for (long i = 0; i < m.size(); i++)
		{
			sum(i) = m[i] / B;
		}
	}
	return sum;
}


inline double& Vector::operator()(size_t i)
{
	return  m[i];
}

inline size_t Vector::size()const
{
	return  m.size();
}


inline Matrix Matrix::Identity(size_t rows, size_t columns) {
	Matrix m(rows, columns);
#pragma omp parallel for
	for (long i = 0; i < rows; i++)
	{
		if (i < columns)
		{
			m(i, i) = 1.0;
		}
	}
	return m;
}
inline Matrix Matrix::Diagonal() {
	size_t u = this->Rows();
	Matrix m1(u, this->Columns());

#pragma omp parallel for
	for (long i = 0; i < u; i++)
	{
		m1(i, i) = (*this)(i, i);
	}
	return m1;
}
inline Vector Matrix::Gaussian_algorithm(Matrix* m, Vector* v)
{
	size_t u = (*m).Rows();
	size_t t = (*m).Columns();
	Matrix m2(u, t + 1);
	Vector result(v->size());

#pragma omp parallel for
	for (long i = 0; i < u; i++)
	{
		for (size_t j = 0; j < t; j++)
		{
			m2(i, j) = (*m)(i, j);
		}
	}

#pragma omp parallel for
	for (long i = 0; i < u; i++)
	{
		m2(i, t) = (*v)(i);
	}

#pragma omp parallel for
	for (long j = 0; j < t; j++)
	{
		for (size_t i = 0; i < t; i++)
		{
			if (i != j)
			{
				double b = m2(i, j) / m2(j, j);
				for (size_t k = 0; k < t + 1; k++)
				{
					m2(i, k) = m2(i, k) - b * m2(j, k);
				}
			}
		}
	}

#pragma omp parallel for
	for (long i = 0; i < t; i++)
	{
		result(i) = m2(i, t) / m2(i, i);
	}

	return result;
}




inline Matrix::Matrix(size_t rows, size_t columns) {
	r = rows;
	c = columns;
	m.resize(r * c,0.0);
}
inline void Matrix::resize(size_t rows, size_t columns)
{
	r = rows;
	c = columns;
	m.resize(r * c, 0.0);
}
inline double& Matrix::operator()(size_t row, size_t column)
{
	return  m[row * c + column];
}


inline Matrix Matrix::operator+(Matrix B) {
	Matrix sum(c, r);
#pragma omp parallel for
	for (long i = 0; i < r; i++)
	{
		for (size_t j = 0; j < c; j++)
		{
			sum(i, j) = (*this)(i, j) + B(i, j);
		}
	}
	return sum;
}

inline Matrix Matrix::operator-(Matrix B) {
	Matrix diff(r, c);
#pragma omp parallel for
	for (long i = 0; i < r; i++)
	{
		for (size_t j = 0; j < c; j++)
		{
			diff(i, j) = (*this)(i, j) - B(i, j);
		}
	}

	return diff;
}

inline Matrix Matrix::operator*(Matrix B) {
	Matrix multip(r, B.Columns());

	if (c == B.Rows())
	{
		#pragma omp parallel for
		for (long i = 0; i < r; i++)
		{
			for (size_t j = 0; j < B.Columns(); j++)
			{
				double temp = 0.0;
				for (size_t k = 0; k < c; k++)
				{
					temp += (*this)(i, k) * B(k, j);
				}
				multip(i, j) = temp;
			}

		}
	}
	return multip;
}

inline Matrix Matrix::operator*(double scalar) {
	Matrix result(r, c);
	#pragma omp parallel for
	for (long i = 0; i < r; i++)
	{
		for (size_t j = 0; j < c; j++)
		{
			result(i, j) = (*this)(i, j) * scalar;
		}
	}
	return result;
}

inline Vector Matrix::operator*(Vector B) {

	Vector result((this)->Rows());
	if (B.size() == (this)->Columns())
	{
		#pragma omp parallel for
		for (long i = 0; i < r; i++)
		{
			double sum = 0;
			for (size_t j = 0; j < c; j++)
			{
				sum += (*this)(i, j) * B(j);
			}
			result(i) = sum;
		}
	}
	return result;
}

inline Matrix Matrix::operator/(double scalar) {
	Matrix result(r, c);
	#pragma omp parallel for
	for (long i = 0; i < r; i++)
	{
		for (size_t j = 0; j < c; j++)
		{
			result(i, j) = (*this)(i, j) / scalar;
		}
	}
	return result;
}

inline size_t Matrix::Rows() const
{
	return this->r;
}

inline size_t Matrix::Columns() const
{
	return this->c;
}

inline Matrix Matrix::Transpose()
{
	Matrix t(c, r);
	#pragma omp parallel for
	for (long i = 0; i < c; i++)
	{
		for (size_t j = 0; j < r; j++) {
			t(i, j) = (*this)(j, i);
		}
	}
	return t;
}

inline void Matrix::printmatrix()
{
	#pragma omp parallel for
	for (long i = 0; i < (*this).Rows(); i++)
	{
		for (size_t j = 0; j < (*this).Columns(); j++)
		{
			cout << (*this)(i, j) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

inline void Vector::printvector()
{
	#pragma omp parallel for
	for (long i = 0; i < m.size(); i++)
	{
		cout << m[i] << " ";

	}
	cout << endl;
}


inline Matrix Jacobian(Vector* X, Vector* beta, bool withbe, bool withtheta)
{
	size_t t = (*X).size();
	Matrix m1(t, 4);
	if (withbe && withtheta)
	{
		m1.resize(t, 5);

		double alpha = (*beta)(0), gamma = (*beta)(1), x0 = (*beta)(2), theta = (*beta)(3);
		#pragma omp parallel for
		for (long i = 0; i < t; i++)
		{
			double h = (*X)(i) - x0;
			double c = pow(h, 2.0);
			double w = alpha * c + gamma;
			double g = pow(w, 2.0);
			double p = g * w;
			double q = exp(-theta / w);
			m1(i, 0) = -((w - theta) * q + w) * c / p;
			m1(i, 1) = ((-w + theta) * q - w) / p;
			m1(i, 2) = (2.0 * ((w - theta) * q + w)) * alpha * h / p;
			m1(i, 3) = -q / g;
			m1(i, 4) = 1.0;
		}
	}
	else if (withtheta &&!withbe)
	{

		double alpha = (*beta)(0), gamma = (*beta)(1), x0 = (*beta)(2), theta = (*beta)(3);
		#pragma omp parallel for
		for (long i = 0; i < t; i++)
		{
			double h = (*X)(i) - x0;
			double c = pow(h, 2.0);
			double w = alpha * c + gamma;
			double g = pow(w, 2.0);
			double p = g * w;
			double q = exp(-theta / w);
			m1(i, 0) = -c*((w - theta) * q + w) / p;
			m1(i, 1) = ((-w + theta) * q - w) / p;
			m1(i, 2) = (2.0 *alpha * h*((w - theta) * q + w)) / p;
			m1(i, 3) = -q / g;
		}
	}
	else if (!withtheta && withbe)
	{
		#pragma omp parallel for
		for (long i = 0; i < t; i++)
		{
			double alpha = (*beta)(0), gamma = (*beta)(1), x0 = (*beta)(2);
			double h = (*X)(i) - x0;
			double c = pow(h, 2.0);
			double w = alpha * c + gamma;
			double g = pow(w, 2.0);
			m1(i, 0) = -c / g;
			m1(i, 1) = -1.0 / g;
			m1(i, 2) = 2.0 * alpha * h / g;
			m1(i, 3) = 1.0;
		}
	}

	return m1;
}

inline Vector f(Vector* X, Vector* beta,bool withbe,bool withtheta)
{
	Vector m((*X).size());
	
	if (withbe && withtheta)
	{
		double alpha = (*beta)(0), gamma = (*beta)(1), x0 = (*beta)(2), theta = (*beta)(3),be= (*beta)(4);
		#pragma omp parallel for
		for (long i = 0; i < (*X).size(); i++)
		{
			double w = alpha * pow((*X)(i) - x0, 2.0) + gamma;
			m(i) = 1.0 / w + exp(-theta / w) / w + be;
		}
	}
	else if(!withbe && withtheta)
	{
		double alpha = (*beta)(0), gamma = (*beta)(1), x0 = (*beta)(2), theta = (*beta)(3);
		#pragma omp parallel for
		for (long i = 0; i < (*X).size(); i++)
		{
			double w = alpha * pow((*X)(i) - x0, 2.0) + gamma;
			m(i) = 1.0 / w + exp(-theta / w) / w;
		}
	}
	else if (withbe && !withtheta)
	{
		double alpha = (*beta)(0), gamma = (*beta)(1), x0 = (*beta)(2), be=(*beta)(3);
		#pragma omp parallel for
		for (long i = 0; i < (*X).size(); i++)
		{
			double w = alpha * pow((*X)(i) - x0, 2.0) + gamma;
			m(i) = 1.0 / w + be;
		}
	}

	return m;
}

inline Vector directionalderivative(Vector* X,Vector*fxbeta, Vector* beta, Vector* delta, Matrix* J,bool withbe,bool withtheta)
{
	double h = 0.1;
	Vector temp = (*beta) + ((*delta) * h);
	return  ((f(X, &temp,withbe,withtheta) - *fxbeta) / h - (*J) * (*delta)) * (2 / h);

}

inline double error(Vector* Y, Vector* X, Vector* beta, bool withbe, bool withtheta, outlier_criterion method, double tolerance, double additionaldata, valarray<bool>* indices)
{
	valarray<double> tmp= (f(X, beta,withbe,withtheta)-(*Y)).array();
	return error(&tmp, method, tolerance, additionaldata, indices);
}
inline double error(Vector* Y, Vector* X, Vector* beta, bool withbe, bool withtheta)
{
	Vector func = f(X, beta, withbe, withtheta)-(*Y);
	return (func*func);
}


inline double error(valarray<double>* errs, outlier_criterion method,double tolerance,double additionaldata,valarray<bool>*indicestouse)
{

	valarray <double> errs1 = *errs;

	

	double m = 0, MAD = 0, average = 0, stdev = 0, S = 0, Q = 0, T = 0, biweightmidvariance = 0, pbmv=0;
	size_t pointnumber = (*errs).size();
	size_t pointnumberhalf = pointnumber / 2;

	//compute the user selected estimator for all the errors
	switch (method)
	{
	case errorfunction_vanishing_tolerance_multiplies_standard_deviation_of_error:
	{
		stdeviation(&errs1, &stdev, &average, pointnumber);
		break;
	}
	case errorfunction_vanishing_tolerance_is_significance_in_Grubbs_test:
	{
		stdeviation(&errs1, &stdev, &average, pointnumber);
		break;
	}
	case errorfunction_vanishing_tolerance_is_decision_in_MAD_ESTIMATION:
	{
		MAD = MAD_estimator(&errs1, &m, pointnumber, pointnumberhalf);
		break;
	}
	case errorfunction_vanishing_tolerance_is_biweight_midvariance:
	{
		biweightmidvariance = onestepbiweightmidvariance(&errs1, &m, pointnumber, pointnumberhalf);
		break;
	}
	
	case errorfunction_vanishing_tolerance_is_decision_in_S_ESTIMATION:
	{
		S = S_estimator(&errs1, pointnumber);
		m = median(&errs1, pointnumber, pointnumberhalf);
		break;
	}
	case errorfunction_vanishing_tolerance_is_decision_in_Q_ESTIMATION:
	{
		m = median(&errs1, pointnumber, pointnumberhalf);
		Q = Q_estimator(&errs1, pointnumber);
		break;
	}
	case errorfunction_vanishing_tolerance_is_decision_in_T_ESTIMATION:
	{
		m = median(&errs1, pointnumber, pointnumberhalf);
		T = T_estimator(&errs1, pointnumber);
		break;
	}
	case errorfunction_vanishing_use_peirce_criterion:
	{
		valarray<double> temp(pointnumber);
		for (size_t w = 0; w < pointnumber; w++)
		{
			temp[w] = errs1[w] * errs1[w];
		}

		stdeviation(&temp, NULL, &average, pointnumber);
		break;
	}


	case errorfunction_linear_tolerance_multiplies_standard_deviation_of_error:
	{
		stdeviation(&errs1, &stdev, &average, pointnumber);
		break;
	}
	case errorfunction_linear_tolerance_is_significance_in_Grubbs_test:
	{
		stdeviation(&errs1, &stdev, &average, pointnumber);
		break;
	}
	case errorfunction_linear_tolerance_is_decision_in_MAD_ESTIMATION:
	{
		MAD = MAD_estimator(&errs1, &m, pointnumber, pointnumberhalf);
		break;
	}
	case errorfunction_linear_tolerance_is_biweight_midvariance:
	{
		biweightmidvariance = onestepbiweightmidvariance(&errs1, &m, pointnumber, pointnumberhalf);
		break;
	}
	
	case errorfunction_linear_tolerance_is_decision_in_S_ESTIMATION:
	{
		S = S_estimator(&errs1, pointnumber);
		m = median(&errs1, pointnumber, pointnumberhalf);
		break;
	}
	case errorfunction_linear_tolerance_is_decision_in_Q_ESTIMATION:
	{
		m = median(&errs1, pointnumber, pointnumberhalf);
		Q = Q_estimator(&errs1, pointnumber);
		break;
	}
	case errorfunction_linear_tolerance_is_decision_in_T_ESTIMATION:
	{
		m = median(&errs1, pointnumber, pointnumberhalf);
		T = T_estimator(&errs1, pointnumber);
		break;
	}
	case errorfunction_linear_use_peirce_criterion:
	{
		valarray<double> temp(pointnumber);
		for (size_t w = 0; w < pointnumber; w++)
		{
			temp[w] = errs1[w] * errs1[w];
		}

		stdeviation(&temp, NULL, &average, pointnumber);
		break;
	}
	}

	double err_entire = 0;
	size_t points = 0;
	for (size_t i = 0; i < errs->size(); i++)
	{
		double err = (*errs)[i];
		switch (method)
		{
		case errorfunction_vanishing_tolerance_is_maximum_squared_error:
		{
			double G = (err) * (err);
			if (G <= fabs(tolerance))
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_vanishing_tolerance_multiplies_standard_deviation_of_error:
		{
			double G = fabs(err - average);
			if (G <= tolerance * stdev)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_vanishing_tolerance_is_significance_in_Grubbs_test:
		{
			double G = fabs(err - average);
			if (G <= additionaldata * stdev)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_vanishing_tolerance_is_decision_in_MAD_ESTIMATION:
		{
			double G = fabs((err - m) / MAD);
			if (G <= tolerance)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_vanishing_tolerance_is_biweight_midvariance:
		{
			double G = fabs((err - m) / biweightmidvariance);
			if (G <= tolerance)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				(*indicestouse)[i] = false;
			}
			break;
		}
	

		case errorfunction_vanishing_tolerance_is_decision_in_S_ESTIMATION:
		{
			double G = (fabs(err - m) / S);
			if (G <= tolerance)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				(*indicestouse)[i] = false;
			}
			break;
		}

		case errorfunction_vanishing_tolerance_is_decision_in_Q_ESTIMATION:
		{
			double G = fabs((err - m) / Q);
			if (G <= tolerance)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_vanishing_tolerance_is_decision_in_T_ESTIMATION:
		{
			double G = fabs((err - m) / T);
			if (G <= tolerance)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_vanishing_use_peirce_criterion:
		{
			double G = (err) * (err);
			if (G <= average * additionaldata)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_linear_tolerance_is_maximum_squared_error:
		{
			double G = (err) * (err);
			if (G <= fabs(tolerance))
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				err_entire += abs(err);
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_linear_tolerance_multiplies_standard_deviation_of_error:
		{
			double G = fabs(err - average);
			if (G <= tolerance * stdev)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				err_entire += abs(err);
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_linear_tolerance_is_significance_in_Grubbs_test:
		{
			double G = fabs(err - average);
			if (G <= additionaldata * stdev)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				err_entire += abs(err);
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_linear_tolerance_is_decision_in_MAD_ESTIMATION:
		{
			double G = fabs((err - m) / MAD);
			if (G <= tolerance)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				err_entire += abs(err);
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_linear_tolerance_is_biweight_midvariance:
		{
			double G = fabs((err - m) / biweightmidvariance);
			if (G <= tolerance)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				err_entire += abs(err);
				(*indicestouse)[i] = false;
			}
			break;
		}
	
		case errorfunction_linear_tolerance_is_decision_in_S_ESTIMATION:
		{
			double G = (fabs(err - m) / S);
			if (G <= tolerance)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
				
			}
			else
			{
				err_entire += abs(err);
				(*indicestouse)[i] = false;
			}
			break;
		}

		case errorfunction_linear_tolerance_is_decision_in_Q_ESTIMATION:
		{
			double G = fabs((err - m) / Q);
			if (G <= tolerance)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				err_entire += abs(err);
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_linear_tolerance_is_decision_in_T_ESTIMATION:
		{
			double G = fabs((err - m) / T);
			if (G <= tolerance)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;				
			}
			else
			{
				err_entire += abs(err);
				(*indicestouse)[i] = false;
			}
			break;
		}
		case errorfunction_linear_use_peirce_criterion:
		{
			double G = (err) * (err);
			if (G <= average * additionaldata)
			{
				err_entire += err * err;
				points++;
				(*indicestouse)[i] = true;
			}
			else
			{
				err_entire += abs(err);
				(*indicestouse)[i] = false;
			}
			break;
		}
		}
	}

	return err_entire/(double)points;
}



inline void computedelta1(Vector* delta, Matrix* G, Vector* v, double lambda)
{
	Matrix tmp = ((*G) + ((*G).Diagonal() * lambda));
	*delta = (*G).Gaussian_algorithm(&tmp, v);
}

inline void computedelta2(Vector* delta, Matrix* G, Vector* v, double lambda)
{
	Matrix tmp = ((*G) + ((*G).Diagonal() * lambda));
	Vector tmp2 = (*v) * (-0.5);
	*delta = (*G).Gaussian_algorithm(&tmp, &tmp2);
}



//computes the factorial
inline double factorial(size_t n)
{
	double ret = 1.00;
	for (size_t i = 2; i <= n; ++i)
		ret *= (double)i;
	return ret;
}
//helper function for student's t distribution 
//from the algorithm in Smiley W. Cheng, James C. Fu, Statistics & Probability Letters 1 (1983), 223-227
inline double H(double y, size_t nu)
{
	double sum = 0;

	for (size_t j = 1; j <= size_t((nu + 1) / 2 - 1); j++)
	{
		sum += factorial((size_t)j) * factorial(size_t(j - 1.0)) / (pow(4, -((double)j)) * factorial((size_t)2.0 * j)) * pow((1 + y * y / nu), -((double)j));
	}
	return y / (2 * sqrt((double)nu)) * sum;
}




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
inline double median(valarray<double> *arr, size_t n, size_t nhalf)
{

#if (__cplusplus == 201703L) && !defined(MACOSX)
	nth_element(std::execution::par, std::begin(*arr), std::begin(*arr) + nhalf, std::begin(*arr) + n);
#else
	nth_element(std::begin(*arr), std::begin(*arr) + nhalf, std::begin(*arr) + n);
#endif


	double  med = (*arr)[nhalf];
	if (n % 2 == 0)
	{
#if __cplusplus == 201703L && !defined(MACOSX)
		auto max_it = max_element(std::execution::par, std::begin(*arr), std::begin(*arr) + nhalf);
#else
		auto max_it = max_element(std::begin(*arr), std::begin(*arr) + nhalf);
#endif

		med = (*max_it + med) / 2.0;
	}
	return med;
}
inline double median(valarray<float>* arr, size_t n, size_t nhalf)
{

#if (__cplusplus == 201703L) && !defined(MACOSX)
	nth_element(std::execution::par, std::begin(*arr), std::begin(*arr) + nhalf, std::begin(*arr) + n);
#else
	nth_element(std::begin(*arr), std::begin(*arr) + nhalf, std::begin(*arr) + n);
#endif


	float  med = (*arr)[nhalf];
	if (n % 2 == 0)
	{
#if __cplusplus == 201703L && !defined(MACOSX)
		auto max_it = max_element(std::execution::par, std::begin(*arr), std::begin(*arr) + nhalf);
#else
		auto max_it = max_element(std::begin(*arr), std::begin(*arr) + nhalf);
#endif

		med = (*max_it + med) / 2.0;
	}
	return med;
}
// computes the lower median of an array. expects the size of the array.
inline double lowmedian(valarray<double> *arr, size_t n)
{
	size_t m = (size_t)(floor(((double)n + 1.0) / 2.0) - 1.0);

#if __cplusplus == 201703L && !defined(MACOSX)
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


inline void hyperbola_linear_regression(double* ferr, valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression, valarray<bool>* usepoint)
{
	*ferr = DBL_MAX;

	valarray<double> x2 = (*x)[*usepoint];
	valarray<double> line_y2 = (*line_y)[*usepoint];
	size_t usedpoints = x2.size();

	valarray<double>line_x2(usedpoints);
	valarray<double> temp(usedpoints);
	valarray<double> temp2(usedpoints);
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
		
			double k;
			for (size_t i = 0; i < usedpoints; i++)
			{
				k = (x2)[i] - (double)h;
				line_x2[i] = k * k;
			}
			for (size_t i = 0; i < usedpoints; i++)
			{
				size_t q = 0;
				for (size_t j = 0; j < usedpoints; j++)
				{
					double t = line_x2[j] - line_x2[i];
					if (t != 0)
					{
						stacks1[q] = (line_y2[j] - line_y2[i]) / t;
						q++;
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

			double thiserr = 0, k2;
			for (size_t n = 0; n < usedpoints; n++)
			{
				k2 =fabs(thisslope * line_x2[n] + thisintercept) - line_y2[n];
				thiserr += k2 * k2;
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
		double yaverage = sumy / (double)usedpoints;
		//for each possible focus position, fit a hyperbola. return the hyperbola with the focus position that has the smallest error to the points 
		for (long h = minfocus; h <= maxfocus; h++)
		{
			double sumx = 0, sumxy = 0, sumxx = 0, k;
			for (size_t n = 0; n < usedpoints; n++)
			{
				k = (x2)[n] - (double)h;
				k = k * k;
				line_x2[n] = k;
				sumx += k;
				sumxy += k * line_y2[n];
				sumxx += k * k;
			}

			double xaverage = sumx / (double)usedpoints;
			double yaverage = sumy / (double)usedpoints;

			double t = sumxx - sumx * xaverage;


			double thisslope = (sumxy - sumx * yaverage) / t;
			double thisintercept = yaverage - thisslope * xaverage;

			double thiserr = 0, k2;


			for (size_t n = 0; n < usedpoints; n++)
			{
				k2 = fabs(thisslope * line_x2[n] + thisintercept) - line_y2[n];
				thiserr += k2 * k2;
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
}

inline void hyperbola_linear_regression(double* ferr, valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression,outlier_criterion method, double tolerance,double additionaldata, valarray<bool>* usepoint)
{
	*ferr = DBL_MAX;


	size_t usedpoints = (*x).size();

	valarray<double>line_x2(usedpoints);
	valarray<double> temp(usedpoints);
	valarray<double> temp2(usedpoints);
	valarray<double> errs(usedpoints);
	valarray<bool> ind(true,usedpoints);
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

			double k;
			for (size_t i = 0; i < usedpoints; i++)
			{
				k = (*x)[i] - (double)h;
				line_x2[i] = k * k;
			}
			for (size_t i = 0; i < usedpoints; i++)
			{
				size_t q = 0;
				for (size_t j = 0; j < usedpoints; j++)
				{
					double t = line_x2[j] - line_x2[i];
					if (t != 0)
					{
						stacks1[q] = ((*line_y)[j] - (*line_y)[i]) / t;
						q++;
					}
				}
				stacks2[i] = median(&stacks1, q, q / 2);
			}

			double thisslope = median(&stacks2, usedpoints, halfsize);

			for (size_t n = 0; n < usedpoints; n++)
			{
				stacki1[n] = (*line_y)[n] - thisslope * line_x2[n];
			}
			double thisintercept = median(&stacki1, usedpoints, halfsize);

			for (size_t n = 0; n < usedpoints; n++)
			{
				errs[n] = fabs(thisslope * line_x2[n] + thisintercept) - (*line_y)[n];
			}

			double thiserr = error(&errs, method, tolerance, additionaldata,&ind);

			if (thiserr < *ferr)
			{
				*fintercept = thisintercept;
				*fslope = thisslope;
				*ferr = thiserr;
				*focpos = h;
				*usepoint = ind;
			}
		}
	}
	else
	{
		//make a linear regression fit for the hyperbola

		double sumy = 0;
		for (size_t n = 0; n < usedpoints; n++)
		{
			sumy += (*line_y)[n];
		}
		double yaverage = sumy / (double)usedpoints;
		//for each possible focus position, fit a hyperbola. return the hyperbola with the focus position that has the smallest error to the points 
		for (long h = minfocus; h <= maxfocus; h++)
		{
			double sumx = 0, sumxy = 0, sumxx = 0, k;
			for (size_t n = 0; n < usedpoints; n++)
			{
				k = (*x)[n] - (double)h;
				k = k * k;
				line_x2[n] = k;
				sumx += k;
				sumxy += k * (*line_y)[n];
				sumxx += k * k;
			}

			double xaverage = sumx / (double)usedpoints;
			double yaverage = sumy / (double)usedpoints;

			double t = sumxx - sumx * xaverage;


			double thisslope = (sumxy - sumx * yaverage) / t;
			double thisintercept = yaverage - thisslope * xaverage;

			for (size_t n = 0; n < usedpoints; n++)
			{
				errs[n] = fabs(thisslope * line_x2[n] + thisintercept) - (*line_y)[n];
			}

			double thiserr = error(&errs, method, tolerance, additionaldata,&ind);

			if (thiserr < *ferr)
			{
				*fintercept = thisintercept;
				*fslope = thisslope;
				*ferr = thiserr;
				*focpos = h;
				*usepoint = ind;
			}
		}
	}


}


inline void powerfunction_linear_regression(double* ferr, valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression,outlier_criterion method, double tolerance,double additionaldata, valarray<bool>* usepoint)
{
	*ferr = DBL_MAX;


	size_t usedpoints = (*x).size();

	valarray<double>line_x2(usedpoints);
	valarray<double> temp(usedpoints);
	valarray<double> temp2(usedpoints);
	valarray<double> errs(usedpoints);
	valarray<bool> ind(true,usedpoints);
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
	
			double k;
			for (size_t i = 0; i < usedpoints; i++)
			{
				k = (*x)[i] - (double)h;
				line_x2[i] = k * k;
			}
			for (size_t i = 0; i < usedpoints; i++)
			{
				size_t q = 0;
				for (size_t j = 0; j < usedpoints; j++)
				{
					double t = line_x2[j] - line_x2[i];
					if (t != 0)
					{
						stacks1[q] = ((*line_y)[j] - (*line_y)[i]) / t;
						q++;
					}
				}
				stacks2[i] = median(&stacks1, q, q / 2);
			}

			double thisslope = median(&stacks2, usedpoints, halfsize);

			for (size_t n = 0; n < usedpoints; n++)
			{
				stacki1[n] = (*line_y)[n] - thisslope * line_x2[n];
			}
			double thisintercept = median(&stacki1, usedpoints, halfsize);

			for (size_t n = 0; n < usedpoints; n++)
			{
				errs[n] = 1.0 / (thisslope * line_x2[n] + thisintercept) - 1.0 / (*line_y)[n];
			}

			double thiserr = error(&errs, method, tolerance, additionaldata,&ind);

			if (thiserr < *ferr)
			{
				*fintercept = thisintercept;
				*fslope = thisslope;
				*ferr = thiserr;
				*focpos = h;
				*usepoint = ind;
			}
		}
	}
	else
	{
		//make a linear regression fit for the hyperbola
		double sumy = 0;
		for (size_t n = 0; n < usedpoints; n++)
		{
			sumy += (*line_y)[n];
		}
		double yaverage = sumy / (double)usedpoints;
		//for each possible focus position, fit a hyperbola. return the hyperbola with the focus position that has the smallest error to the points 
		for (long h = minfocus; h <= maxfocus; h++)
		{
			double sumx = 0, sumxy = 0, sumxx = 0,k;
			for (size_t n = 0; n < usedpoints; n++)
			{
				k = (*x)[n] - (double)h;
				k = k * k;
				line_x2[n] = k;
				sumx += k;
				sumxy += k * (*line_y)[n];
				sumxx += k * k;
			}

			double xaverage = sumx / (double)usedpoints;
			double yaverage = sumy / (double)usedpoints;

			double t = sumxx - sumx * xaverage;
			

			double thisslope = (sumxy - sumx * yaverage) / t;
			double thisintercept = yaverage - thisslope * xaverage;


			for (size_t n = 0; n < usedpoints; n++)
			{
				errs[n] = 1.0 / (thisslope * line_x2[n] + thisintercept) - 1.0 / (*line_y)[n];
			}

			double thiserr = error(&errs, method, tolerance, additionaldata,&ind);

			if (thiserr < *ferr)
			{
				*fintercept = thisintercept;
				*fslope = thisslope;
				*ferr = thiserr;
				*focpos = h;
				*usepoint = ind;
			}
		}
	}

}


inline void powerfunction_linear_regression(double* ferr, valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, long* focpos, double* fintercept, double* fslope, bool use_median_regression, valarray<bool>* usepoint)
{
	*ferr = DBL_MAX;

	valarray<double> x2 = (*x)[*usepoint];
	valarray<double> line_y2 = (*line_y)[*usepoint];

	size_t usedpoints = x2.size();

	valarray<double>line_x2(usedpoints);
	valarray<double> temp(usedpoints);
	valarray<double> temp2(usedpoints);
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

			double k;
			for (size_t i = 0; i < usedpoints; i++)
			{
				k = (x2)[i] - (double)h;
				line_x2[i] = k * k;
			}
			for (size_t i = 0; i < usedpoints; i++)
			{
				size_t q = 0;
				for (size_t j = 0; j < usedpoints; j++)
				{
					double t = line_x2[j] - line_x2[i];
					if (t != 0)
					{
						stacks1[q] = (line_y2[j] - line_y2[i]) / t;
						q++;
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

			double thiserr = 0, k2;
			for (size_t n = 0; n < usedpoints; n++)
			{
				k2 = 1.0 / (thisslope * line_x2[n] + thisintercept) - 1.0 / line_y2[n];
				thiserr += k2 * k2;
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
		double yaverage = sumy / (double)usedpoints;
		//for each possible focus position, fit a hyperbola. return the hyperbola with the focus position that has the smallest error to the points 
		for (long h = minfocus; h <= maxfocus; h++)
		{
			double sumx = 0, sumxy = 0, sumxx = 0, k;
			for (size_t n = 0; n < usedpoints; n++)
			{
				k = (x2)[n] - (double)h;
				k = k * k;
				line_x2[n] = k;
				sumx += k;
				sumxy += k * line_y2[n];
				sumxx += k * k;
			}

			double xaverage = sumx / (double)usedpoints;
			double yaverage = sumy / (double)usedpoints;

			double t = sumxx - sumx * xaverage;


			double thisslope = (sumxy - sumx * yaverage) / t;
			double thisintercept = yaverage - thisslope * xaverage;

			double thiserr = 0, k2;

			for (size_t n = 0; n < usedpoints; n++)
			{
				k2 = 1.0 / (thisslope * line_x2[n] + thisintercept) - 1.0 / line_y2[n];
				thiserr += k2 * k2;
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
}


//starts the regression function and ensures that if the minimum point of the hyperbola is close to minfocus or maxfocus, the regression is started again on an enlarged area at the side close to minfocus or maxfocus. The size of this area
//is determined by the scale parameter.
inline void hyperbola_regression_extendsearcharea(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression, valarray<bool>* indicestouse)
{

	hyperbola_linear_regression(err, x, line_y, minfocus, maxfocus, focpos, fintercept, fslope, use_median_regression, indicestouse);
	

	long dist = abs(*focpos - maxfocus);
	if ((dist <= 10) && (scale > 1))
	{

		const long middle = lround(((double)maxfocus + (double)minfocus) / 2.0);
		minfocus = maxfocus;
		maxfocus = (long)((double)middle + ((double)maxfocus - (double)middle) * fabs(scale));
		long thisfocpos;
		double thiserr, thisslope, thisintercept;
		hyperbola_linear_regression(&thiserr, x, line_y, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope, use_median_regression, indicestouse);
		

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
		hyperbola_linear_regression(&thiserr, x, line_y, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope, use_median_regression, indicestouse);
		
		if ((thiserr < *err))
		{
			*err = thiserr;
			*fintercept = thisintercept;
			*fslope = thisslope;
			*focpos = thisfocpos;
		}
	}
}

inline void hyperbola_regression_extendsearcharea(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression, outlier_criterion method, double tolerance,double additionaldata, valarray<bool>* indicestouse)
{

	hyperbola_linear_regression(err, x, line_y, minfocus, maxfocus, focpos, fintercept, fslope, use_median_regression,method, tolerance, additionaldata, indicestouse);


	long dist = abs(*focpos - maxfocus);
	if ((dist <= 10) && (scale > 1))
	{

		const long middle = lround(((double)maxfocus + (double)minfocus) / 2.0);
		minfocus = maxfocus;
		maxfocus = (long)((double)middle + ((double)maxfocus - (double)middle) * fabs(scale));
		long thisfocpos;
		double thiserr, thisslope, thisintercept;
		hyperbola_linear_regression(&thiserr, x, line_y, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope, use_median_regression,  method, tolerance, additionaldata, indicestouse);


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
		hyperbola_linear_regression(&thiserr, x, line_y, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope, use_median_regression, method, tolerance, additionaldata, indicestouse);

		if ((thiserr < *err))
		{
			*err = thiserr;
			*fintercept = thisintercept;
			*fslope = thisslope;
			*focpos = thisfocpos;
		}
	}
}

//starts the regression function and ensures that if the minimum point of the hyperbola is close to minfocus or maxfocus, the regression is started again on an enlarged area at the side close to minfocus or maxfocus. The size of this area
//is determined by the scale parameter.
inline void powerfunction_regression_extendsearcharea(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression, valarray<bool>* indicestouse)
{

	powerfunction_linear_regression(err, x, line_y, minfocus, maxfocus, focpos,fintercept, fslope, use_median_regression, indicestouse);


	long dist = abs(*focpos - maxfocus);
	if ((dist <= 10) && (scale > 1))
	{

		const long middle = lround(((double)maxfocus + (double)minfocus) / 2.0);
		minfocus = maxfocus;
		maxfocus = (long)((double)middle + ((double)maxfocus - (double)middle) * fabs(scale));
		long thisfocpos;
		double thiserr, thisslope, thisintercept;
		powerfunction_linear_regression(&thiserr, x, line_y, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope,use_median_regression, indicestouse);


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
		powerfunction_linear_regression(&thiserr, x, line_y, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope,use_median_regression, indicestouse);

		if ((thiserr < *err))
		{
			*err = thiserr;
			*fintercept = thisintercept;
			*fslope = thisslope;
			*focpos = thisfocpos;
		}
	}
}

inline void powerfunction_regression_extendsearcharea(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, double* err, double* fintercept, double* fslope, long* focpos, bool use_median_regression, outlier_criterion method, double tolerance,double additionaldata, valarray<bool>* indicestouse)
{

	powerfunction_linear_regression(err, x, line_y, minfocus, maxfocus, focpos, fintercept, fslope, use_median_regression, method, tolerance,additionaldata, indicestouse);


	long dist = abs(*focpos - maxfocus);
	if ((dist <= 10) && (scale > 1))
	{

		const long middle = lround(((double)maxfocus + (double)minfocus) / 2.0);
		minfocus = maxfocus;
		maxfocus = (long)((double)middle + ((double)maxfocus - (double)middle) * fabs(scale));
		long thisfocpos;
		double thiserr, thisslope, thisintercept;
		powerfunction_linear_regression(&thiserr, x, line_y, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope, use_median_regression,method, tolerance, additionaldata, indicestouse);


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
		powerfunction_linear_regression(&thiserr, x, line_y, minfocus, maxfocus, &thisfocpos, &thisintercept, &thisslope, use_median_regression,method, tolerance, additionaldata, indicestouse);

		if ((thiserr < *err))
		{
			*err = thiserr;
			*fintercept = thisintercept;
			*fslope = thisslope;
			*focpos = thisfocpos;
		}
	}
}




//helper function for student's t distribution from the algorithm in
// Smiley W. Cheng, James C. Fu, Statistics & Probability Letters 1 (1983), 223-227
inline double G(double y, size_t nu)
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
	double diff1 = pointnumber - numberofoutliers;

	double res = 0.0;
	
	double Q = pow(numberofoutliers, (numberofoutliers / pointnumber)) * pow(diff1, (diff1 / pointnumber)) / pointnumber;
	double a = 1.0, b = 0.0; 
	while (abs(a - b) > (pointnumber * 2.0e-15)) 
	{
		double l = pow(a, numberofoutliers);
		if (l == 0.0)
		{
			l = 1.0e-5;
		}		
		res = 1.0 + (diff1 - fittingparameters) / numberofoutliers * (1.0 - pow(pow((pow(Q, pointnumber) / l), (1.0 / diff1)), 2.0));
		if (res >= 0)
		{
			b = a;
			a = exp((res - 1.0) / 2.0) * erfc(sqrt(res) / sqrt(2.0));
		}
		else
		{
			b = a;
			res = 0.0;
		}
	}
	return res;
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

#if __cplusplus == 201703L && !defined(MACOSX)
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
#if __cplusplus == 201703L && !defined(MACOSX)
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

inline void hyperbola_findmodel(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, valarray<bool>* usedpoint, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method, size_t pointnumber, size_t pointnumberhalf)
{
	long this_focpos = 0;
	double thiserr = DBL_MAX, thisslope = 0.0, thisintercept = 0.0;
	//makes a fit with the initial usedpoint array as a mask
	hyperbola_regression_extendsearcharea(x, line_y, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &this_focpos, use_median_regression, usedpoint);



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
			maybe_inliner o = {p, z	};
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
		case Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error:
		{
			stdeviation(&err, &stdev, &average, pointnumber);
			break;
		}
		case Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test:
		{
			stdeviation(&err, &stdev, &average, pointnumber);
			break;
		}
		case Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION:
		{
			MAD = MAD_estimator(&err, &m, pointnumber, pointnumberhalf);
			break;
		}
		case Least_trimmed_squares_tolerance_is_biweight_midvariance:
		{
			biweightmidvariance = onestepbiweightmidvariance(&err, &m, pointnumber, pointnumberhalf);
			break;
		}
		
		case Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION:
		{
			S = S_estimator(&err, pointnumber);
			m = median(&err, pointnumber, pointnumberhalf);
			break;
		}
		case Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION:
		{
			m = median(&err, pointnumber, pointnumberhalf);
			Q = Q_estimator(&err, pointnumber);
			break;
		}
		case Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION:
		{
			m = median(&err, pointnumber, pointnumberhalf);
			T = T_estimator(&err, pointnumber);
			break;
		}
		case Least_trimmed_squares_use_peirce_criterion:
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
			case Least_trimmed_squares_tolerance_is_maximum_squared_error:
			{
				double G = mp[j].error * mp[j].error;
				if (G > fabs(tolerance))
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error:
			{
				double G = fabs(mp[j].error - average);
				if (G > tolerance * stdev)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test:
			{
				double G = fabs(mp[j].error - average);
				if (G > additionaldata * stdev)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / MAD);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_biweight_midvariance:
			{
				double G = fabs((mp[j].error - m) / biweightmidvariance);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			

			case Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION:
			{
				double G = (fabs(mp[j].error - m) / S);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}

			case Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / Q);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / T);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_use_peirce_criterion:
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
}

inline void powerfunction_findmodel(valarray<double>* x, valarray<double>* line_y, long minfocus, long maxfocus, double scale, valarray<bool>* usedpoint, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method, size_t pointnumber, size_t pointnumberhalf)
{
	long this_focpos = 0;
	double thiserr = DBL_MAX, thisslope = 0.0, thisintercept = 0.0,thistheta=0.0;
	//makes a fit with the initial usedpoint array as a mask
	powerfunction_regression_extendsearcharea(x, line_y, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &this_focpos,use_median_regression, usedpoint);



	vector<maybe_inliner> mp;
	mp.reserve(pointnumber);

	valarray<double>err(pointnumber);

	//makes a loop over all points
	for (size_t p = 0; p < pointnumber; p++)
	{

		double xh = ((double)(*x)[p] - (double)this_focpos);
		xh *= xh;
		//measures the error between the initial fit and point p
		double z = 1.0/(thisslope * xh + thisintercept) - 1.0/(*line_y)[p];

		//store the error and the squared error
		err[p] = z;

		//if point p is removed by the mask, add it as a maybe inlier
		if (!(*usedpoint)[p])
		{
			maybe_inliner o={p,z};
			//o.point = p;
			//o.error = z;
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
		case Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error:
		{
			stdeviation(&err, &stdev, &average, pointnumber);
			break;
		}
		case Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test:
		{
			stdeviation(&err, &stdev, &average, pointnumber);
			break;
		}
		case Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION:
		{
			MAD = MAD_estimator(&err, &m, pointnumber, pointnumberhalf);
			break;
		}
		case Least_trimmed_squares_tolerance_is_biweight_midvariance:
		{
			biweightmidvariance = onestepbiweightmidvariance(&err, &m, pointnumber, pointnumberhalf);
			break;
		}
		
		case Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION:
		{
			S = S_estimator(&err, pointnumber);
			m = median(&err, pointnumber, pointnumberhalf);
			break;
		}
		case Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION:
		{
			m = median(&err, pointnumber, pointnumberhalf);
			Q = Q_estimator(&err, pointnumber);
			break;
		}
		case Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION:
		{
			m = median(&err, pointnumber, pointnumberhalf);
			T = T_estimator(&err, pointnumber);
			break;
		}
		case Least_trimmed_squares_use_peirce_criterion:
		{
			valarray<double> temp(pointnumber);
			for (size_t w = 0; w < pointnumber; w++)
			{
				temp[w] = err[w] * err[w];
			}

			stdeviation(&temp, NULL, &average, pointnumber);
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
			case Least_trimmed_squares_tolerance_is_maximum_squared_error:
			{
				double G = mp[j].error * mp[j].error;
				if (G > fabs(tolerance))
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error:
			{
				double G = fabs(mp[j].error - average);
				if (G > tolerance * stdev)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test:
			{
				double G = fabs(mp[j].error - average);
				if (G > additionaldata * stdev)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / MAD);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_biweight_midvariance:
			{
				double G = fabs((mp[j].error* mp[j].error - m) / biweightmidvariance);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
		

			case Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION:
			{
				double G = (fabs(mp[j].error - m) / S);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}

			case Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / Q);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / T);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_use_peirce_criterion:
			{
				double G = err[j]*err[j];
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
}

inline void powerfunction_nonlinear_regression(double* err, valarray<double>* x, valarray<double>* y, double* alpha, double* gamma, double* b, double* theta, long* focpos,outlier_criterion method, double tolerance, double additionaldata, valarray<bool>* usedpoint, long minfocus, long maxfocus, double scale, size_t stop_after_iterations, double stop_after_seconds, bool withbe, bool withtheta)

{

	Vector q((*x).size());
	Vector Y((*y).size());
	Vector X((*x).size());
	valarray<bool> usedpoints1(true,x->size());
	valarray<bool> usedpoints2(true,x->size());
	for (size_t w = 0; w < X.size(); w++)
	{
		if ((*usedpoint)[w])
		{
			X(w) = (*x)[w];
			Y(w) = 1.0 / (*y)[w];
		}
	}
	Vector beta(4);
	Vector delta1(4);
	Vector delta2(4);
	Vector delta1a(4);
	Vector delta2a(4);
	Vector delta1b(4);
	Vector delta2b(4);

	if (withbe && withtheta)
	{
		beta.resize(5);
		delta1.resize(5);
		delta2.resize(5);
		delta1a.resize(5);
		delta2a.resize(5);
		delta1b.resize(5);
		delta2b.resize(5);

		beta(0) = *alpha;
		beta(1) = *gamma;
		beta(2) = *focpos;
		beta(3) = *theta;
		beta(4) = *b;
	}
	else if (withtheta && !withbe)
	{
		beta(0) = *alpha;
		beta(1) = *gamma;
		beta(2) = *focpos;
		beta(3) = *theta;
	}
	else if (withbe && !withtheta)
	{
		beta(0) = *alpha;
		beta(1) = *gamma;
		beta(2) = *focpos;
		beta(3) = *b;
	}

	*err = DBL_MAX;



	double seconds = 0.0;
	double lambda = 4.0;
	double increment = 1.5;
	double decrement = 5.0;


	int count = 0;
	const long middle = lround(((double)maxfocus + (double)minfocus) / 2.0);
	long maxfocus1 = maxfocus;
	long minfocus1 = minfocus;
	if (scale <= 1.0)
	{
		scale = 1.0;
	}

	maxfocus = (long)((double)middle + ((double)maxfocus1 - (double)middle) * fabs(scale));
	minfocus = (long)((double)middle - ((double)middle - (double)minfocus1) * fabs(scale));

	auto start = std::chrono::steady_clock::now();

	do
	{
		Matrix J = Jacobian(&X, &beta, withbe, withtheta);


		Matrix Jt = J.Transpose();

		Vector fxbeta = f(&X, &beta, withbe, withtheta);
		Vector  v = Jt * (Y - fxbeta);

		Matrix G = Jt * J;

		computedelta1(&delta1a, &G, &v, lambda);

		double  norm1a = delta1a * delta1a;
		if (isnan(norm1a))
		{
			break;
		}


		delta1 = delta1a;


		Vector v1b = Jt * directionalderivative(&X, &fxbeta, &beta, &delta1a, &J, withbe, withtheta);

		computedelta2(&delta1b, &G, &v1b, lambda);

		double  norm1b = delta1b * delta1b;
		if ((2.0 * sqrt(norm1b) / sqrt(norm1a)) < 0.1)
		{
			delta1 = delta1 + delta1b;
		}
		norm1a = delta1a * delta1a;


		Vector temp = beta + delta1;

		double  s1 = error(&Y, &X, &temp, withbe, withtheta,method, tolerance,additionaldata,&usedpoints1);


		computedelta1(&delta2a, &G, &v, lambda / decrement);

		double  norm2a = delta2a * delta2a;

		if (isnan(norm2a))
		{
			break;
		}

		delta2 = delta2a;

		Vector v2b = Jt * directionalderivative(&X, &fxbeta, &beta, &delta2a, &J, withbe, withtheta);

		computedelta2(&delta2b, &G, &v2b, lambda / decrement);

		double  norm2b = delta2b * delta2b;
		if ((2.0 * sqrt(norm2b) / sqrt(norm2a)) < 0.1)
		{
			delta2 = delta2 + delta2b;
		}

		norm2a = delta2a * delta2a;


		Vector temp2 = beta + delta2;


		double  s2 = error(&Y, &X, &temp2, withbe, withtheta,method,tolerance,additionaldata, &usedpoints2);


		if (norm2a < 0.01)
		{
			break;
		}
		if (norm1a < 0.01)
		{
			break;
		}
	

		if ((s1 > *err) && (s2 > *err))
		{
			lambda *= increment;
			continue;
		}

		if (s2 < s1)
		{
			lambda = lambda / decrement;
			beta = temp2;
			*usedpoint = usedpoints2;
			*err = s2;
			if (fabs(delta2(2)) < 0.01)
			{
				break;
			}
		}
		else
		{
			beta = temp;
			*usedpoint = usedpoints1;
			*err = s1;
			if (fabs(delta1(2)) < 0.01)
			{
				break;
			}
		}

		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		seconds = elapsed_seconds.count();

		count++;
		if (count > stop_after_iterations)
			break;

	} while (seconds < stop_after_seconds);

	*alpha = beta(0);
	*gamma = beta(1);
	*focpos =(long) round(beta(2));

	if (withtheta)
		*theta = beta(3);

	if (withbe && !withtheta)
		*b = beta(3);

	if (withbe && withtheta)
		*b = beta(4);

	if (*focpos > maxfocus)
	{
		*focpos = maxfocus;
	}
	if (*focpos < minfocus)
	{
		*focpos = minfocus;
	}
}



inline void powerfunction_nonlinear_regression(double* err, valarray<double>* x, valarray<double>* y, double* alpha, double* gamma,double *b, double* theta, long* focpos, valarray<bool>* usedpoint, long minfocus, long maxfocus, double scale, size_t stop_after_iterations, double stop_after_seconds, bool withbe, bool withtheta)

{

	
	valarray<double> x2 = (*x)[*usedpoint];
	valarray<double> y2 = (*y)[*usedpoint];

	Vector q((x2).size());
	Vector Y((y2).size());
	Vector X((x2).size());

	for (size_t w = 0; w < X.size(); w++)
	{
		if((*usedpoint)[w])
		{
			X(w) = x2[w];
			Y(w) = 1.0/y2[w];
		}
	}
	Vector beta(4);
	Vector delta1(4);
	Vector delta2(4);
	Vector delta1a(4);
	Vector delta2a(4);
	Vector delta1b(4);
	Vector delta2b(4);

	if (withbe && withtheta)
	{
		beta.resize(5);
		delta1.resize(5);
		delta2.resize(5);
		delta1a.resize(5);
		delta2a.resize(5);
		delta1b.resize(5);
		delta2b.resize(5);

		beta(0) = *alpha;
		beta(1) = *gamma;
		beta(2) = *focpos;
		beta(3) = *theta;
		beta(4) = *b;
	}
	else if (withtheta && !withbe)
	{
		beta(0) = *alpha;
		beta(1) = *gamma;
		beta(2) = *focpos;
		beta(3) = *theta;
	}
	else if  (withbe && !withtheta)
	{
		beta(0) = *alpha;
		beta(1) = *gamma;
		beta(2) = *focpos;
		beta(3) = *b;
	}

	*err = DBL_MAX;



	double seconds = 0.0;
	double lambda = 4.0;
	double increment = 1.5;
	double decrement = 5.0;

	
	int count = 0;
	const long middle = lround(((double)maxfocus + (double)minfocus) / 2.0);
	long maxfocus1 = maxfocus;
	long minfocus1 = minfocus;
	if (scale <= 1.0)
	{
		scale = 1.0;
	}

	maxfocus = (long)((double)middle + ((double)maxfocus1 - (double)middle) * fabs(scale));
	minfocus = (long)((double)middle - ((double)middle - (double)minfocus1) * fabs(scale));

	auto start = std::chrono::steady_clock::now();

	do
	{
		Matrix J = Jacobian(&X, &beta,withbe,withtheta);


		Matrix Jt = J.Transpose();

		Vector fxbeta = f(&X, &beta,withbe,withtheta);
		Vector  v = Jt * (Y -fxbeta);

		Matrix G = Jt * J;

		computedelta1(&delta1a, &G, &v, lambda);

		double  norm1a = delta1a * delta1a;
		if (isnan(norm1a))
		{
			break;
		}


		delta1 = delta1a;


		Vector v1b = Jt * directionalderivative(&X,&fxbeta, &beta, &delta1a, &J,withbe,withtheta);

		computedelta2(&delta1b, &G, &v1b, lambda);

		double  norm1b = delta1b * delta1b;
		if ((2 * sqrt(norm1b) / sqrt(norm1a)) < 0.1)
		{
			delta1 = delta1 + delta1b;
		}
		norm1a = delta1a * delta1a;


		Vector temp(beta.size());
		temp = beta + delta1;

		double  s1 = error(&Y, &X, &temp,withbe,withtheta);


		computedelta1(&delta2a, &G, &v, lambda / decrement);

		double  norm2a = delta2a * delta2a;

		if (isnan(norm2a))
		{
			break;
		}

		delta2 = delta2a;

		Vector v2b = Jt * directionalderivative(&X,&fxbeta, &beta, &delta2a, &J,withbe,withtheta);

		computedelta2(&delta2b, &G, &v2b, lambda / decrement);

		double  norm2b = delta2b * delta2b;
		if ((2 * sqrt(norm2b) / sqrt(norm2a)) < 0.1)
		{
			delta2 = delta2 + delta2b;
		}

		norm2a = delta2a * delta2a;

		Vector temp2(beta.size());
		temp2 = beta + delta2;

	
		double  s2 = error(&Y, &X, &temp2,withbe,withtheta);
		

		if (norm2a < 0.01)
		{
			break;
		}
		if (norm1a < 0.01)
		{
			break;
		}

		
		if ((s1 > *err) && (s2 > *err))
		{
			lambda *= increment;
			continue;
		}

		if (s2 < s1)
		{		
			lambda = lambda / decrement;
			beta = temp2;
			
			*err = s2;
		}
		else
		{
			beta = temp;
			*err = s1;
		}

		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		seconds = elapsed_seconds.count();

		count++;
		if (count > stop_after_iterations)
			break;

	} while (seconds<stop_after_seconds);
	
	 *alpha = beta(0);
	 *gamma = beta(1);
	 *focpos= (long) beta(2);

	 if (withtheta)
	 {
		 *theta = beta(3);
	 }

	 if (withbe && !withtheta)
	 {
		 *b = beta(3);
	 }
	 if (withbe && withtheta)
	 {
		 *b = beta(4);
	 }


	 *err = *err / x2.size();


	 if (*focpos > maxfocus)
	 {
		 *focpos = maxfocus;
	 }
	 if (*focpos < minfocus)
	 {
		*focpos = minfocus;
	 }
}


inline void powerfunction_findmodel(valarray<double>* x, valarray<double>* line_y, double alpha, double gamma, long focpos, long minfocus, long maxfocus, double scale, valarray<bool>* usedpoint, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method, size_t pointnumber, size_t pointnumberhalf, size_t stop_after_iterations, double stop_after_seconds,bool withbe, bool withtheta)

{
	double thiserr = DBL_MAX;
	double beta = 0,theta=10*gamma;

	powerfunction_nonlinear_regression(&thiserr,x, line_y, &alpha, &gamma, &beta, &theta, &focpos, usedpoint,minfocus,maxfocus,scale,  stop_after_iterations,  stop_after_seconds,withbe, withtheta);



	vector<maybe_inliner> mp;
	mp.reserve(pointnumber);

	valarray<double>err(pointnumber);
	valarray<double>err_sqr(pointnumber);
	//makes a loop over all points
	for (size_t p = 0; p < pointnumber; p++)
	{

		double xh = ((double)(*x)[p] - (double)focpos);
		xh *= xh;
		double c = (alpha * xh + gamma);
		double z = 0;
		if (withbe && withtheta)
		{
			z = 1.0 / c + exp(-theta / c) / c + beta - 1.0 / (*line_y)[p];
		}
		else if (!withbe && withtheta)
		{
			 z = 1.0 / c + exp(-theta / c) / c - 1.0 / (*line_y)[p];
		}
		else if (withbe && !withtheta)
		{
			z = 1.0 / c + beta - 1.0 / (*line_y)[p];
		}

		//store the error and the squared error
		err[p] = z;
		err_sqr[p] = z * z;
		//if point p is removed by the mask, add it as a maybe inlier
		if (!(*usedpoint)[p])
		{
			maybe_inliner o={p,z};
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
		case Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error:
		{
			stdeviation(&err, &stdev, &average, pointnumber);
			break;
		}
		case Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test:
		{
			stdeviation(&err, &stdev, &average, pointnumber);
			break;
		}
		case Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION:
		{
			MAD = MAD_estimator(&err, &m, pointnumber, pointnumberhalf);
			break;
		}
		case Least_trimmed_squares_tolerance_is_biweight_midvariance:
		{
			biweightmidvariance = onestepbiweightmidvariance(&err, &m, pointnumber, pointnumberhalf);
			break;
		}
		
		case Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION:
		{
			S = S_estimator(&err, pointnumber);
			m = median(&err, pointnumber, pointnumberhalf);
			break;
		}
		case Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION:
		{
			m = median(&err, pointnumber, pointnumberhalf);
			Q = Q_estimator(&err, pointnumber);
			break;
		}
		case Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION:
		{
			m = median(&err, pointnumber, pointnumberhalf);
			T = T_estimator(&err, pointnumber);
			break;
		}
		case Least_trimmed_squares_use_peirce_criterion:
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
			case Least_trimmed_squares_tolerance_is_maximum_squared_error:
			{
				double G = mp[j].error * mp[j].error;
				if (G > fabs(tolerance))
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error:
			{
				double G = fabs(mp[j].error - average);
				if (G > tolerance * stdev)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test:
			{
				double G = fabs(mp[j].error - average);
				if (G > additionaldata * stdev)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / MAD);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_biweight_midvariance:
			{
				double G = fabs((mp[j].error - m) / biweightmidvariance);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			

			case Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION:
			{
				double G = (fabs(mp[j].error - m) / S);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}

			case Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / Q);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION:
			{
				double G = fabs((mp[j].error - m) / T);
				if (G > tolerance)
				{
					isoutlier = true;
				}
				break;
			}
			case Least_trimmed_squares_use_peirce_criterion:
			{
				double G = err[j] * err[j];
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
}

//the function which is exported by the library. It uses the ransac to compute a hyperbolic fit where it selects points as inliers and outliers. The documentation is provided in the header file.
inline void hyperbola_Least_trimmed_squares(double *error,long*focpos,double *slope, double *intercept, valarray<bool> *returnindices,size_t pointnumber,size_t pointnumberhalf, size_t maximum_number_of_outliers, size_t stop_after_numberofiterations_without_improvement,double stop_after_seconds,valarray<bool> *indices, valarray<double>* xv, valarray<double> *line_yv,long minfocus, long maxfocus, double scale, double tolerance, double additionaldata,bool use_median_regression, outlier_criterion rejection_method)
{

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
			arr[i] = *indices;
			std::next_permutation(std::begin(*indices), std::end(*indices));
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
#if __cplusplus == 201703L && !defined(MACOSX)

		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				hyperbola_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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
			});
		std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				hyperbola_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
				mtx.lock();
				if (thiserr < *error)
				{
					*error = thiserr;
					*focpos = thisfocpos;
					*slope = thisslope;
					*intercept = thisintercept;
					*returnindices = arri;
				}
				mtx.unlock();
			});
#else
#pragma omp parallel for
		for (long i = 0; i < (long)numbercomp; i++)
		{
			hyperbola_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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

#pragma omp parallel for
		for (long i = 0; i < (long)count2; i++)
		{
			double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
			long thisfocpos = 0;
			hyperbola_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
#pragma omp critical
			{
				if (thiserr <*error)
				{
					*error = thiserr;
					*focpos = thisfocpos;
					*slope = thisslope;
					*intercept = thisintercept;
					*returnindices = arr2[i];
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
				arr[i] = *indices;
				std::next_permutation(std::begin(*indices), std::end(*indices));
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
#if __cplusplus == 201703L && !defined(MACOSX)
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					hyperbola_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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
				});

			std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
				{
					double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
					long thisfocpos = 0;
					hyperbola_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
					mtx.lock();
					if (thiserr < *error)
					{
						*error = thiserr;
						*focpos = thisfocpos;
						*slope = thisslope;
						*intercept = thisintercept;
						*returnindices = arri;
					}
					mtx.unlock();
				});
#else			
#pragma omp parallel for
			for (long i = 0; i < (long)number_of_attempts; i++)
			{
				hyperbola_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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

#pragma omp parallel for
			for (long i = 0; i < count2; i++)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				hyperbola_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
#pragma omp critical
				{
					if (thiserr < *error)
					{
						*error = thiserr;
						*focpos = thisfocpos;
						*slope = thisslope;
						*intercept = thisintercept;
						*returnindices = arr2[i];
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
			arr[i] = *indices;
			std::next_permutation(std::begin(*indices), std::end(*indices));
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
#if __cplusplus == 201703L && !defined(MACOSX)
		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				hyperbola_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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
			});

		std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				hyperbola_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
				mtx.lock();
				if (thiserr < *error)
				{
					*error = thiserr;
					*focpos = thisfocpos;
					*slope = thisslope;
					*intercept = thisintercept;
					*returnindices = arri;
				}
				mtx.unlock();
			});
#else
#pragma omp parallel for
		for (long i = 0; i < (long)s; i++)
		{
			hyperbola_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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

#pragma omp parallel for
		for (long i = 0; i < count2; i++)
		{
			double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
			long thisfocpos = 0;
			hyperbola_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
#pragma omp critical
			{
				if (thiserr < *error)
				{
					*error = thiserr;
					*focpos = thisfocpos;
					*slope = thisslope;
					*intercept = thisintercept;
					*returnindices = arr2[i];
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
				shuffle(std::begin(*indices), std::end(*indices), urng);
				arr[i] = *indices;
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
#if __cplusplus == 201703L && !defined(MACOSX)
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					hyperbola_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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
						hyperbola_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
						mtx.lock();
						if (thiserr < *error)
						{
							*error = thiserr;
							*focpos = thisfocpos;
							*slope = thisslope;
							*intercept = thisintercept;
							*returnindices = arri;
							//if we had found a better model set counter 1 =0
							counter1 = 0;
						}
						mtx.unlock();
					}
				});
#else
#pragma omp parallel for
			for (long i = 0; i < number_of_attempts; i++)
			{
				hyperbola_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
				//for each model we tried to find increase counter1 by 1
#pragma omp atomic
				counter1++;
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

#pragma omp parallel for
			for (long i = 0; i < count2; i++)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				hyperbola_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
#pragma omp critical
				{
					if (thiserr < *error)
					{
						*error = thiserr;
						*focpos = thisfocpos;
						*slope = thisslope;
						*intercept = thisintercept;
						*returnindices = arr2[i];
						//if we had found a better model set the counter 1 to zero
						counter1 = 0;
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


}

//the function which is exported by the library. It uses the ransac to compute a hyperbolic fit where it selects points as inliers and outliers. The documentation is provided in the header file.
inline void powerfunction_Least_trimmed_squares(double* error, long* focpos, double* slope, double* intercept, valarray<bool>* returnindices, size_t pointnumber, size_t pointnumberhalf, size_t maximum_number_of_outliers, size_t stop_after_numberofiterations_without_improvement, double stop_after_seconds, valarray<bool>* indices, valarray<double>* xv, valarray<double>* line_yv, long minfocus, long maxfocus, double scale, double tolerance, double additionaldata, bool use_median_regression, outlier_criterion rejection_method)
{

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
			arr[i] = *indices;
			std::next_permutation(std::begin(*indices), std::end(*indices));
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
#if __cplusplus == 201703L && !defined(MACOSX)

		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				powerfunction_findmodel(xv, line_yv,minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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
			});
		std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0;
				long thisfocpos = 0;
				powerfunction_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept,&thisslope,&thisfocpos, use_median_regression, &arri);
				mtx.lock();
				if (thiserr < *error)
				{
					*error = thiserr;
					*focpos = thisfocpos;
					*slope = thisslope;
					*intercept = thisintercept;
					*returnindices = arri;
				}
				mtx.unlock();
			});
#else
#pragma omp parallel for
		for (long i = 0; i < (long)numbercomp; i++)
		{
			powerfunction_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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

#pragma omp parallel for
		for (long i = 0; i < (long)count2; i++)
		{
			double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0,thistheta=0;
			long thisfocpos = 0;
			powerfunction_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
#pragma omp critical
			{
				if (thiserr < *error)
				{
					*error = thiserr;
					*focpos = thisfocpos;
					*slope = thisslope;
					*intercept = thisintercept;
					*returnindices = arr2[i];
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
				arr[i] = *indices;
				std::next_permutation(std::begin(*indices), std::end(*indices));
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
#if __cplusplus == 201703L && !defined(MACOSX)
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					powerfunction_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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
				});

			std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
				{
					double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0,thistheta=0;
					long thisfocpos = 0;
					powerfunction_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
					mtx.lock();
					if (thiserr < *error)
					{
						*error = thiserr;
						*focpos = thisfocpos;
						*slope = thisslope;
						*intercept = thisintercept;
						*returnindices = arri;
					}
					mtx.unlock();
				});
#else			
#pragma omp parallel for
			for (long i = 0; i < (long)number_of_attempts; i++)
			{
				powerfunction_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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

#pragma omp parallel for
			for (long i = 0; i < count2; i++)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0,thistheta=0;
				long thisfocpos = 0;
				powerfunction_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
#pragma omp critical
				{
					if (thiserr < *error)
					{
						*error = thiserr;
						*focpos = thisfocpos;
						*slope = thisslope;
						*intercept = thisintercept;
						*returnindices = arr2[i];
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
			arr[i] = *indices;
			std::next_permutation(std::begin(*indices), std::end(*indices));
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
#if __cplusplus == 201703L && !defined(MACOSX)
		std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
			{
				powerfunction_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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
			});

		std::for_each(std::execution::par, std::begin(arr2), std::begin(arr2) + count2, [&](valarray<bool>& arri)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0, thistheta = 0;
				long thisfocpos = 0;
				powerfunction_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
				mtx.lock();
				if (thiserr < *error)
				{
					*error = thiserr;
					*focpos = thisfocpos;
					*slope = thisslope;
					*intercept = thisintercept;
					*returnindices = arri;
				}
				mtx.unlock();
			});
#else
#pragma omp parallel for
		for (long i = 0; i < (long)s; i++)
		{
			powerfunction_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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

#pragma omp parallel for
		for (long i = 0; i < count2; i++)
		{
			double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0,thistheta=0;
			long thisfocpos = 0;
			powerfunction_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
#pragma omp critical
			{
				if (thiserr < *error)
				{
					*error = thiserr;
					*focpos = thisfocpos;
					*slope = thisslope;
					*intercept = thisintercept;
					*returnindices = arr2[i];
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
				shuffle(std::begin(*indices), std::end(*indices), urng);
				arr[i] = *indices;
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
#if __cplusplus == 201703L && !defined(MACOSX)
			std::for_each(std::execution::par, std::begin(arr), std::end(arr), [&](valarray<bool>& arri)
				{
					powerfunction_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arri, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
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
						double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0,thistheta=0;
						long thisfocpos = 0;
						powerfunction_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arri);
						mtx.lock();
						if (thiserr < *error)
						{
							*error = thiserr;
							*focpos = thisfocpos;
							*slope = thisslope;
							*intercept = thisintercept;
							*returnindices = arri;
							//if we had found a better model set counter 1 =0
							counter1 = 0;
						}
						mtx.unlock();
					}
				});
#else
#pragma omp parallel for
			for (long i = 0; i < number_of_attempts; i++)
			{
				powerfunction_findmodel(xv, line_yv, minfocus, maxfocus, scale, &arr[i], tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf);
				//for each model we tried to find increase counter1 by 1
#pragma omp atomic
				counter1++;
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

#pragma omp parallel for
			for (long i = 0; i < count2; i++)
			{
				double thiserr = DBL_MAX, thisslope = 0, thisintercept = 0,thistheta=0;
				long thisfocpos = 0;
				powerfunction_regression_extendsearcharea(xv, line_yv, minfocus, maxfocus, scale, &thiserr, &thisintercept, &thisslope, &thisfocpos, use_median_regression, &arr2[i]);
#pragma omp critical
				{
					if (thiserr < *error)
					{
						*error = thiserr;
						*focpos = thisfocpos;
						*slope = thisslope;
						*intercept = thisintercept;
						*returnindices = arr2[i];
						//if we had found a better model set the counter 1 to zero
						counter1 = 0;
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


}


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
	if (rejection_method == Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test)
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
	valarray <double> xv(pointnumber);
	for (size_t i = 0; i < pointnumber; i++)
	{
		if (i >= maximum_number_of_outliers)
		{
			indices[i] = true;
		}

		line_yv[i] = y[i]*y[i];
		xv[i] = (double)(x[i]);
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
	if (rejection_method == Least_trimmed_squares_use_peirce_criterion)
	{

		additionaldata = peirce(pointnumber, maximum_number_of_outliers, 3);
	}
	//generate additional data for the Grubbs test
	if (rejection_method == Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test)
	{
		additionaldata = crit(tolerance, pointnumber);
	}
	if (rejection_method == errorfunction_vanishing_use_peirce_criterion)
	{
		additionaldata = peirce(pointnumber, maximum_number_of_outliers, 3);
	}
	//generate additional data for the Grubbs test
	if (rejection_method == errorfunction_vanishing_tolerance_is_significance_in_Grubbs_test)
	{
		additionaldata = crit(tolerance, pointnumber);
	}
	if (rejection_method == errorfunction_linear_use_peirce_criterion)
	{
		additionaldata = peirce(pointnumber, maximum_number_of_outliers, 3);
	}
	//generate additional data for the Grubbs test
	if (rejection_method == errorfunction_linear_tolerance_is_significance_in_Grubbs_test)
	{
		additionaldata = crit(tolerance, pointnumber);
	}


	if (rejection_method < 10)
	{
		if (maximum_number_of_outliers > 0)
		{
			hyperbola_Least_trimmed_squares(&error, focpos, &slope, &intercept, &indices2, pointnumber, pointnumberhalf, maximum_number_of_outliers, stop_after_numberofiterations_without_improvement, stop_after_seconds, &indices, &xv, &line_yv, minfocus, maxfocus, scale, tolerance, additionaldata, use_median_regression, rejection_method);
		}
		else
		{
			hyperbola_regression_extendsearcharea(&xv, &line_yv, minfocus, maxfocus, scale, &error, &intercept, &slope, focpos, use_median_regression, &indices);
			indices2 = indices;
		}
	}
	else
	{
		hyperbola_regression_extendsearcharea(& xv, &line_yv, minfocus, maxfocus, scale, & error, & intercept, & slope, focpos, use_median_regression, rejection_method, tolerance, additionaldata, & indices2);
		indices2 = indices;
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
//returns the inverse power of an image in fourier mode.
double image::invpower()
{
	return pinvpower;
}
//returns the power of an image in fourier mode.
double image::power()
{
	return ppower;
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
		ppower=datafilling0(fptr, focuser_position);
		pinvpower = 1 / ppower;
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
	ppower = datafilling0(fptr, focuser_position);
	pinvpower = 1 / ppower;
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
			return fouriertransform(NULL, &pimagedata_double, dimension1, dimension2,bitpix);
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
			return fouriertransform(&pimagedata_float, NULL, dimension1, dimension2,bitpix);
		}
	}

}

//analyzes the fourier data of the image and calculates the power. The image data can be either in fload or double.
double image::fouriertransform(vector<float>* p4, vector<double>* p5, size_t dimension1, size_t dimension2,int bitpix)
{
	double outburnedconst;
	if (bitpix == 8)
		outburnedconst = 255;
	else if (bitpix == 16)
		outburnedconst = 65535;
	else if (bitpix == 32)
		outburnedconst = 4294967295;
	else if (bitpix == 64)
		outburnedconst = 18446744073709551615;
	else if (bitpix == -32)
		outburnedconst = 3.402823466e38;
	else if (bitpix == -64)
		outburnedconst = 1.7976931348623158e308;

	if (p4 != NULL)
	{
		//compute the dft of the image
		Mat I( dimension2, dimension1, CV_32F, (*p4).data());

		double minVal;
		double maxVal;
		Point minLoc;
		Point maxLoc;
		
		minMaxLoc(I, &minVal, &maxVal, &minLoc, &maxLoc);
		if ((maxVal == outburnedconst) && (bitpix!=8))
		{
			double c = outburnedconst * 0.01 + minVal;
			for (size_t i = 0; i < dimension2; i++)
			{
				for (size_t j = 0; j < dimension1; j++)
				{
					if (I.at<float>(i, j) >= c)
						I.at<float>(i, j) = c/2;
				}
			}
		}


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

		return power[0] + power[1] + power[2] + power[3];
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

		double minVal;
		double maxVal;
		Point minLoc;
		Point maxLoc;

		minMaxLoc(I, &minVal, &maxVal, &minLoc, &maxLoc);

		if ((maxVal == outburnedconst) && (bitpix != 8))
		{
			double c = outburnedconst * 0.01 + minVal;
			for (size_t i = 0; i < dimension2; i++)
			{
				for (size_t j = 0; j < dimension1; j++)
				{
					if (I.at<double>(i, j) >= c)
						I.at<double>(i, j) = c / 2.0;
				}
			}
		}
		
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



		return power[0] + power[1] + power[2] + power[3];
	}

}

//constructs the image class from a double array.Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t width, size_t height, long focuser_position, vector<double> * imagedata, double hfd, double fwhm)
{
	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (width <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("width was <=2");
	}
	if (height <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("height was <=2");
	}
	if (imagedata->size() ==0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((width * height) != imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}
	ppower=fouriertransform(NULL, imagedata, width, height,-64);
	pinvpower = 1 / ppower;
}
//constructs the image class from a float array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t width, size_t height, long focuser_position, vector<float> * imagedata, double hfd, double fwhm)
{
	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (width <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("width was <=2");
	}
	if (height <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("height was <=2");
	}
	if (imagedata->size() == 0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((width * height)!= imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}
	ppower = fouriertransform(imagedata, NULL, width, height,-32);
	pinvpower = 1 / ppower;
}
//constructs the image class from a long long array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t width, size_t height, long focuser_position, vector<long long> * imagedata, double hfd, double fwhm)
{
	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (width <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("width was <=2");
	}
	if (height <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("height was <=2");
	}
	if (imagedata->size() == 0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((width * height) != imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}

	vector<float>floatvec(imagedata->begin(), imagedata->end());
	ppower = fouriertransform(&floatvec, NULL, width, height,64);
	pinvpower = 1 / ppower;
}
//constructs the image class from a long array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t width, size_t height, long focuser_position, vector<long> * imagedata, double hfd, double fwhm)
{

	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (width <=2)
	{
		pstatus = -1;
		throw std::invalid_argument("width was <=2");
	}
	if (height <=2)
	{
		pstatus = -1;
		throw std::invalid_argument("height was <=2");
	}
	if (imagedata->size() == 0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((width * height) != imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}


	vector<float>floatvec(imagedata->begin(), imagedata->end());
	ppower = fouriertransform(&floatvec, NULL, width, height,32);
	pinvpower = 1 / ppower;
}
//constructs the image class from a short array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t width, size_t height, long focuser_position, vector<short> * imagedata, double hfd, double fwhm)
{
	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (width <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("width was <=2");
	}
	if (height <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("height was <=2");
	}
	if (imagedata->size() == 0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((width * height) != imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}

	
	vector<float>floatvec(imagedata->begin(), imagedata->end());
	ppower = fouriertransform(&floatvec, NULL, width, height,16);
	pinvpower = 1 / ppower;
}
//constructs the image class from an int8_t array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
image::image(size_t width, size_t height, long focuser_position, vector<int8_t> * imagedata, double hfd, double fwhm)
{

	pstatus = 0;
	pfocuser_position = focuser_position;
	phfd = hfd;
	pfwhm = fwhm;
	if (width <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("width was <=2");
	}
	if (height <= 2)
	{
		pstatus = -1;
		throw std::invalid_argument("height was <=2");
	}
	if (imagedata->size() == 0)
	{
		pstatus = -1;
		throw std::invalid_argument("no image data");
	}
	if ((width * height) != imagedata->size())
	{
		pstatus = -1;
		throw std::invalid_argument("dimensions do not correspond to image data size");
	}


	vector<float>floatvec(imagedata->begin(), imagedata->end());
	ppower = fouriertransform(&floatvec, NULL, width, height,8);
	pinvpower = 1 / ppower;

}

//computes the focus point from a vector of validly constructed image classes which must all have status 0.
//the other parameters are as in focusposition_Regression with 2 exceptions: double* main_slope, double* main_intercept are the slope 
// and intercept for the functions invpower=slope(x-focus_point)^2+intercept, where invpower is the given by the invpower method 
//if x is the focus motor position of an image class. It is the inverse of the power function from the fourier analysis.
//also, the parameters  from focusposition_Regression where the fitted curve is returned in a coordinate system where it is represented 
//by a line are omitted. 

bool focusposition_Regression2(std::vector<image>* images, long* focpos, double* main_error, double* alpha, double* beta, double* gamma, double* theta,
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
	if (rejection_method == Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test)
	{
		if (tolerance >= 0.9)
		{
			rejection_method = no_rejection;
		}
	}

	long minfocus = LONG_MAX, maxfocus = LONG_MIN;


	valarray<bool> indices(false, pointnumber);

	valarray<bool> indices2(true, pointnumber);
	valarray<double>yv(pointnumber);
	valarray <double> xv(pointnumber);
	//for all image classes which were successfully created, set an initial bitmask, set the y values to the inverse power, find the minimum and maximum focus value
	double maxbeta = DBL_MAX;
	for (size_t i = 0; i < pointnumber; i++)
	{
		if ((*images)[i].status() == 0)
		{

			if (i >= maximum_number_of_outliers)
			{
				indices[i] = true;
			}

			yv[i] = (*images)[i].invpower();
			xv[i] = (double)(*images)[i].focuser_position();
			if (xv[i] < minfocus)
			{
				minfocus = (long)xv[i];
			}
			if (xv[i] > maxfocus)
			{
				maxfocus = (long)xv[i];
			}
		}
	}

	//check if garbage data was supplied
	if (maxfocus - minfocus == 0)
		return false;


	double  error = DBL_MAX, pgamma = 0, palpha = 0, additionaldata = 0;


	//generate additional data for the Peirce criterion
	if (rejection_method == Least_trimmed_squares_use_peirce_criterion)
	{
		size_t e = 3;
		if (beta != NULL) e++;
		if (theta != NULL) e++;
		additionaldata = peirce(pointnumber, maximum_number_of_outliers, e);
	}
	//generate additional data for the Grubbs test
	if (rejection_method == Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test)
	{
		additionaldata = crit(tolerance, pointnumber);
	}
	if (rejection_method == errorfunction_vanishing_use_peirce_criterion)
	{
		size_t e = 3;
		if (beta != NULL) e++;
		if (theta != NULL) e++;
		additionaldata = peirce(pointnumber, maximum_number_of_outliers, e);
	}
	//generate additional data for the Grubbs test
	if (rejection_method == errorfunction_vanishing_tolerance_is_significance_in_Grubbs_test)
	{
		additionaldata = crit(tolerance, pointnumber);
	}
	if (rejection_method == errorfunction_linear_use_peirce_criterion)
	{
		size_t e = 3;
		if (beta != NULL) e++;
		if (theta != NULL) e++;
		additionaldata = peirce(pointnumber, maximum_number_of_outliers, e);
	}
	//generate additional data for the Grubbs test
	if (rejection_method == errorfunction_linear_tolerance_is_significance_in_Grubbs_test)
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

	bool withbe = false, withtheta = false;

	if (beta != NULL)
	{
		withbe = true;
	}
	if (theta != NULL)
	{
		withtheta = true;
	}

	double ptheta;
	double pbeta;


	if (rejection_method < 10)
	{
		if (maximum_number_of_outliers > 0)
		{
			powerfunction_Least_trimmed_squares(&error, focpos, &palpha, &pgamma, &indices2, pointnumber, pointnumberhalf, maximum_number_of_outliers, stop_after_numberofiterations_without_improvement, stop_after_seconds, &indices, &xv, &yv, minfocus, maxfocus, scale, tolerance, additionaldata, use_median_regression, rejection_method);
		}
		else
		{
			powerfunction_regression_extendsearcharea(&xv, &yv, minfocus, maxfocus, scale, &error, &pgamma, &palpha, focpos, use_median_regression, &indices);
			indices2 = indices;
		}

		if (withbe || withtheta)
		{
			bool hasoutliers = false;
			for (size_t l = 0; l < indices2.size(); l++)
			{
				if (indices2[l] == false)
				{
					hasoutliers = true;
					break;
				}
			}
			if (hasoutliers)
			{
				powerfunction_findmodel(&xv, &yv, palpha, pgamma, *focpos, minfocus, maxfocus, scale, &indices2, tolerance, additionaldata, use_median_regression, rejection_method, pointnumber, pointnumberhalf, stop_after_numberofiterations_without_improvement / 3, stop_after_seconds / 3, withbe, withtheta);

			}
			ptheta = 10.0 * pgamma;
			pbeta = 0.0;
			powerfunction_nonlinear_regression(&error, &xv, &yv, &palpha, &pgamma, &pbeta, &ptheta, focpos, &indices2, minfocus, maxfocus, scale, stop_after_numberofiterations_without_improvement / 3, stop_after_seconds / 3, withbe, withtheta);
		}
	}
	else
	{
		powerfunction_regression_extendsearcharea(&xv, &yv, minfocus, maxfocus, scale, &error, &pgamma, &palpha, focpos, use_median_regression, rejection_method,tolerance,additionaldata, &indices);
		indices2 = indices;

		if (withbe || withtheta)
		{
			ptheta = 10.0 * pgamma;
			pbeta = 0.0;
			powerfunction_nonlinear_regression(&error, &xv, &yv, &palpha, &pgamma, &pbeta, &ptheta, focpos, rejection_method, tolerance, additionaldata, &indices2, minfocus, maxfocus, scale, stop_after_numberofiterations_without_improvement / 3, stop_after_seconds / 3, withbe, withtheta);
		}
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

	if (alpha != NULL)
	{
		*alpha = palpha ;

	}

	if (gamma != NULL)
	{
		*gamma = pgamma ;
	}
	if (theta != NULL)
	{
		*theta = ptheta;
	}

	if (beta != NULL)
	{
		*beta = pbeta;
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
	vector<image>* images1, vector<image>* images2, double* main_error1 , double* alpha1 , double* beta1 , double* gamma1 , double* theta1, vector<size_t>* indicesofusedpoints1 , vector<size_t>* indicesofremovedpoints1 ,
	double* main_error2 , double* alpha2 , double* beta2 , double* gamma2 , double* theta2 , vector<size_t>* indicesofusedpoints2 , vector<size_t>* indicesofremovedpoints2 ,
	double stop_after_seconds , size_t stop_after_numberofiterations_without_improvement, double scale , bool use_median_regression,
	size_t maximum_number_of_outliers , outlier_criterion rejection_method , double tolerance )

{

	if (backslash == NULL)
	{
		return false;
	}

	long focpos1 = 0, focpos2 = 0;



	//make the first fit;
	if (!focusposition_Regression2(images1, &focpos1, main_error1, alpha1, beta1,theta1, gamma1, indicesofusedpoints1, indicesofremovedpoints1, stop_after_seconds,
		stop_after_numberofiterations_without_improvement, 0, scale, use_median_regression, maximum_number_of_outliers, rejection_method, tolerance))
	{
		return false;
	}
	//make the second fit;
	if (!focusposition_Regression2(images2, &focpos2, main_error2, alpha2, beta2, theta2, gamma2, indicesofusedpoints2, indicesofremovedpoints2, stop_after_seconds,
		stop_after_numberofiterations_without_improvement, 0, scale, use_median_regression, maximum_number_of_outliers, rejection_method, tolerance))
	{
		return false;
	}


	//return the computed backslash
	*backslash = focpos2 - focpos1;
	return true;
}
