
// This library implements the algorithms and procedures described inMonthly Notices of the Royal Astronomical Society,
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
// implementation of a Levenberg-Marquardt algorithm in order to fit the fourier data to higher precision
// implementation of robust curve fitting options where the the curve fit is made with all points but the loss function is either 
// linearized to (abs(err)) or set to zero for outliers.



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
// responsible for data taking and various tests of the fourier analysis algorithm


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

// In order to increase the precision of the fit of fourier data, the library has the ability to use a Levenberg-Marquardt algorithm.
// Extensive information about the Levenberg-Marquardt algorithm and some techniques used by this application in order to improve convergence can be found in
// Transtrum, Mark K; Sethna, James P(2012). "Improvements to the Levenberg-Marquardt algorithm for nonlinear least-squares minimization".arXiv:1201.5885




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

// This ia the header file that should be used by external programs which call the library. 


#pragma once
#include "focusinterpolation_exports.h"


using namespace std;

#include <float.h>
#include <fitsio.h>
#include <string> 
#include <vector>


typedef enum
{
	no_rejection,
	Least_trimmed_squares_tolerance_is_maximum_squared_error,
	Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error,
	Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test,
	Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION,
	Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION,
	Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION,
	Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION,
	Least_trimmed_squares_use_peirce_criterion,
	Least_trimmed_squares_tolerance_is_biweight_midvariance,
	errorfunction_vanishing_tolerance_is_maximum_squared_error,
	errorfunction_vanishing_tolerance_multiplies_standard_deviation_of_error,
	errorfunction_vanishing_tolerance_is_significance_in_Grubbs_test,
	errorfunction_vanishing_tolerance_is_decision_in_MAD_ESTIMATION,
	errorfunction_vanishing_tolerance_is_decision_in_S_ESTIMATION,
	errorfunction_vanishing_tolerance_is_decision_in_Q_ESTIMATION,
	errorfunction_vanishing_tolerance_is_decision_in_T_ESTIMATION,
	errorfunction_vanishing_use_peirce_criterion,
	errorfunction_vanishing_tolerance_is_biweight_midvariance,
	errorfunction_linear_tolerance_is_maximum_squared_error,
	errorfunction_linear_tolerance_multiplies_standard_deviation_of_error,
	errorfunction_linear_tolerance_is_significance_in_Grubbs_test,
	errorfunction_linear_tolerance_is_decision_in_MAD_ESTIMATION,
	errorfunction_linear_tolerance_is_decision_in_S_ESTIMATION,
	errorfunction_linear_tolerance_is_decision_in_Q_ESTIMATION,
	errorfunction_linear_tolerance_is_decision_in_T_ESTIMATION,
	errorfunction_linear_use_peirce_criterion,
	errorfunction_linear_tolerance_is_biweight_midvariance,
	errorfunction_linear_tolerance_is_percentagebased_midvariance
}outlier_criterion;


// The function focusposition_Regression is only there for compatibility reasons.

// By now, it is replaced by the function focusposition_Regression2. The image class contains code for image analysis  in fourier space. From a vector of these
// classes, focusposition_Regression2 can interpolate the motor value of the focus point of a telescope.

// In contrast to this new method, focusposition_Regression needs hfd values from another image analysis to interpolates 
// the optimal focuser position. However, both focusposition_Regression2 and focusposition_Regression have very similar parameters and 
// they function in an essentially similar way.

// So we begin the introduction of these functions with the old focusposition_Regression.

// focusposition_Regression interpolates the focus point from a fit with symmetric hyperbolas. It has the ability to solve the trimmed least squares problem
// in order to remove outliers. If specified, it can also use repeaded median regression.
// The algorithm for the linear regression was first published by Stephen King (username STEVE333),
// at https://aptforum.com/phpbb/viewtopic.php?p=25998#p25998

// focusposition_Regression expects the following arguments:

// a vector <long> x with motor positions at which the hfd number was measured. 
// a vector <double> y. The hfd value of a star or a star field corresponding to motor position x [ i ] is assumed to be y [ i ] .
// both x and y should have at least 4 points.

// If the function terminates successfully, the pointer focpos contains the estimated focusposition as a variable of type long that the function computes. This pointer must not be NULL if the
// function is to terminate successfully

// main_error is a pointer to a double value. If the pointer is not NULL, the value of main_error will be the sum of the squares of the differences between a square of the the best fit and a square of the measurement data. 
// Yes, the errors involve two squares because we calculate the squared error of a fit where the hyperbola was converted to a line.


// main_slope and main_intercept are pointers to double values. If the pointers are not NULL, their values will be the slope and intercept of a line given by
// Y=slope*(x-focpos)^2+intercept. The square-root f=sqrt(Y) is the fitted hyperbola.

// indices_of_used_points is a pointer to a vector which, if not NULL, will contain the indices i of the points in x [i] and y [i] that were used for the fit.

// usedpoints_line_x and usedpoints_line_y are pointers to two vectors for the datatype double. If the pointers are not NULL, the vectors contain the points which are used in the fit
// in a coordinate system where the hyperbola is a line. One can plot them together with the line  given by Y=slope*(x-focpos)^2+intercept

// indices_of_removedpoints is a pointer to a vector which, if not NULL, will contain the indices i of the points in x [i] and y [i] that were not used for the fit

// removedpoints_line_x and removedpoints_line_y are pointers to two vectors. If the pointers are not NULL, the vectors contain the points which were not used in the fit
// in the coordinate system where the hyperbola is a line. One can plot them together with the line  given by Y=slope*(x-focpos)^2+intercept



// stop_after_seconds is a parameter that stops the algorithm after a given time in seconds has elapsed.
// stop_after_numberofiterations_without_improvement is a parameter that lets the algorithm stop after it has iterated by stop_after_numberofiterations_without_improvement iterations
// without a further improvement of the error. Note that this parameter is not the iteration number, but it is the number of iterations without further improvement.
// 
// The parameters stop_after_seconds and stop_after_numberofiterations_without_improvement are only used if the binomial coefficient n choose k is larger than 100*(22 choose 11) == 70543200.


// backslash is a parameter that can contain the focuser backslash in steps. The best focus position is corrected with respect to this backslash. If you already have taken account of
// the focuser backslash, for example by setting a suitable overshoot or a final_inwards_movement in APT or a different software or hardware correction of the backslash, set this parameter to 0


// scale is a parameter of the type double that specifies the size of the interval of motor positions where the best focusposition is searched for. 
// The default of scale is 1.
// 
// Sometimes, the focus point may be outside of the interval of motor positions where the hfd was measured.
// let  middle =(max + min) / 2 and max and min be the maximum and minimum motorposition where a hfd was measured.
// If an initial search finds the best focus point within 10 positions to be at the right edge of the measurement interval of motor positions,
// then, if scale>1, the right side of the interval where the best focus is searched is enlarged. The new right boundary is then given by
// max = middle + (max - middle) * abs(scale)
// Similarly, if an initial search finds the best focus point within 10 positions to be on the left side of the measurement interval, the search interval is enlarged, with the new left boundary given by
// min = (middle - (middle - min) * abs(scale).

// use_median_regression is a parameter that specifies whether the RANSAC uses a simple linear regression or a median regression.
// Repeated median regression is slightly more stable against small outliers if one does not use the RANSAC algorithm.


// rejection_method  is a parameter that specifies the method which is used to reject outliers.
// 
// The library can solve the least trimmed squares problem. If such an estimator is selected, the application behaves as follows:
// 
// maximum_number_of_outliers is a parameter that specifies how many outliers the ransac can maximally throw away.

// The algorithm works by first selecting combination of points and fitting them to a hyperbola. 

// 
// Assume you have n datapoints. The algorithm works by searching through either all or (if the binomial coefficient of points over the number of outliers is larger than 20 over 10) randomly generated so - called minimal combinations
// of m=n - maximum_number_of_outliers points.
// 
// An initial fit for a hyperbola is made with these m selected points.  
// 
// Afterwards, this initial hyperbola is then corrected by contributions from other points. A point outside a combination m is added
// if its error from the initial fit is deemed not to be an outlier based on various statistical methods. 
// 
// A new fit with the added point is then made and the process is repeated with another initial combination of points until the best combination had been found.
// The algorithm works in parallel with several threads. So it benefits from processors with multiple cores.
// // 
// 
// The error between the fit w of a minimal combination and a measurement at a motor position x is given by err_p=p(x)-w(x) where w is the squared hfd at x. 

// If

// rejection_method==no_rejection, then the function uses every point for the fit.


// rejection_method==Least_trimmed_squares_tolerance_is_maximum_squared_error, 

// then a point p outside of the  minimal combination are only added to the final fit if its squared error  err_p*err_p fulfills

// err_p*err_p<=abs(tolerance).

// The algorithm  then computes an overall fit with this combination, and then constructs a new set based on a different minimal mode, until the best combination of points was found. 

// To specify the tolerable error directly is useful if the largest tolerable error is known, e.g. from a measurement of the seeing. By taking a series of images, 
// an application may measure the random deviation of the hfd from the average that stars have. With outlier detection methods, one can remove the outliers of large amplitudes.
// setting tolerance_in_sigma_units=false, one can supply a maximally tolerable deviation from the average hfd directly to the library as the tolerance parameter.
// The tolerance value then corresponds to the absolute value of the maximally tolerable hfd deviation.
// note that the error is given with respect to the squared hfd's. I.e if you measure a hfd u at a motor position x, then square its value to g(x)=u(x)^2 and compute the average w(x) of these values. 
// The maximum error is then given by the maximum value of M=abs(w(x)-g(x)), where g(x) is such that it maximizes M. The maximum squared error is given by M^2.

// If 
// rejection_method=Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error,

// then a point p outside of the best minimal combination is added to the final fit if  error err_p fulfills:

// (err_p-a)<= tolerance*s. 

// where a and s are the average and standard deviation of the errors err_p for each measured point with respect to the fit of the given minimal combination w that the algorithm tries.

// If most of your data is good, then a is rather small and you may want to include most points. 
// In that case, set tolerance=1, for example. This excludes only the largest errors that arose from seeing.

// However, assume that for some motor positions far away from the focus point, the hfd values do not follow a hyperbola at all and therefore deviate much. 
// Or assume that you have a large seeing error and you need to throw many points away.

// In that case, your data has many points which have a very large error or difference when compared to a hyperbola. And only a subset of the data close to the focus point may have a small error. 
// In order to exclude larger datasets which do not correspond to a hyperbola at all, you would only want to retain errors which are smaller than the average error a. 
// In that case, set tolerance=-1 or even -2 if many of your points do have very large errors.


// The average and standard deviation is not robust against outliers. If rejection_method== tolerance_is_standard_deviation_of_squared_error_from_mean,
// is used, the tolerance value would have to change depending on whether your data contains many outliers or not.

// Therefore, the library provides more robust methods for accurate removal of outliers where such changes have to happen less frequently.

// If 

// rejection_method==Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION, or 
// rejection_method==Least_trimmed_squares_tolerance_is_biweight_midvariance, or
// rejection_method==Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION, or  
// rejection_method==Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION, or
// rejection_method==Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION, 
// rejection_method==Least_trimmed_squares_tolerance_is_percentagebased_midvariance, 

// then, MAD, biweight_midvariance, S, Q, T or percentage based midvariance estimators are used. 

// A point is then added to a minimal combination if

// abs(err_p-median(err_p))/ estimator <= tolerance

// where median(err_p) is the median of the errors, and estimator is then the MAD, Q, S or T estimator. 
// 
// The MAD,S,Q and T estimators are extensively described in
// Peter J.Rousseeuw, Christophe Croux, Alternatives to the Median - Absolute Deviation
// J. of the Amer. Statistical Assoc. (Theory and Methods), 88 (1993),p. 1273,
// Q and S estimators are better suited for asymmetric distributions than the MAD estimator and the Q estimator is better optimized for small sample sizes than the S estimator.

// The library also can make use of the biweight midvariance estimator that was described in 
// T. C. Beers, K. Flynn and K. Gebhardt,  Astron. J. 100 (1),32 (1990)

// if MAD, S, Q, T estimators, biweight midvariance estimators are used, the tolerance parameter should be around 2...3.


// if 
// rejection_method==Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test,
//
// then a point p is added to a minimal combination if its error err_p is not an outlier according to a modified Grubb's test.
// The tolerance parameter is here a significance value based on student's distribution.
// This means, in order to remove an outlier with 90% significance with Grubb's method, one should set tolerance=0.1, for 80% significance, one should select 0.2.

// In the original version of Grubb's test, one searches the value err_p where |err_p-a|, with a as average is at its maximum.
// One then makes Grubb's test with this value (which depends on the number of samples), and if it is an outlier, one removes this value and repeats the procedure.
// 
// This library modifies Grubb's test a bit, in that it does not include points which fail Grubb's test. It does not search for the largest outlier, removes it 
// and then recalculates the test statistic with an updated sample size. This is done for computational speed, and it should correspond to the behavior of Grubb's test with a very large sample size
// 
// if 
// rejection_method==Least_trimmed_squares_use_peirce_criterion, 

// then a point is added to a minimal combination if it does fulfill Peirce's criterion for not being an outlier. 

// For this the rejection method according to Peirce, the tolerance parameter is not used. 
// Note that number of outliers that Peirce's criterion finds is often strongly influenced by the parameter maximum_number_of_outliers

// tolerance is a parameter that is used for the different rejection methods.
 
// The trimmed least squares problem is computationally intensive. It increases with the binomial coefficient of the number of points over the maximum number of outliers.
// Therefore, if more than 24 points are used, the least trimmed squares prolem is often computationally too expensive.

// The library has a different option. It can perform a linear regression with all points and minimize the squared error, but one
// can choose to linearize errors from points which are deemed to be outliers. This is done with the estimators 

// errorfunction_linear_tolerance_is_maximum_squared_error,
// errorfunction_linear_tolerance_multiplies_standard_deviation_of_error,
// errorfunction_linear_tolerance_is_significance_in_Grubbs_test,
// errorfunction_linear_tolerance_is_decision_in_MAD_ESTIMATION,
// errorfunction_linear_tolerance_is_decision_in_S_ESTIMATION,
// errorfunction_linear_tolerance_is_decision_in_Q_ESTIMATION,
// errorfunction_linear_tolerance_is_decision_in_T_ESTIMATION,
// errorfunction_linear_use_peirce_criterion,
// errorfunction_linear_tolerance_is_biweight_midvariance,
// errorfunction_linear_tolerance_is_percentagebased_midvariance

// The estimators are the same as before in the least trimmed squares approach. 
// This computation is nevertheless different from the trimmed least squares regression, since all points are still
// used for the computation of the slope and intercept parameters
// and just the error function that is used to find the minimum of the hyperbola considers the contributions of
// as linear influences. 

// similarly, one can chose to set errors that are deemed to be outliers to zero. This is done by the estimators
// 
// errorfunction_linear_tolerance_is_maximum_squared_error,
// errorfunction_linear_tolerance_multiplies_standard_deviation_of_error,
// errorfunction_linear_tolerance_is_significance_in_Grubbs_test,
// errorfunction_linear_tolerance_is_decision_in_MAD_ESTIMATION,
// errorfunction_linear_tolerance_is_decision_in_S_ESTIMATION,
// errorfunction_linear_tolerance_is_decision_in_Q_ESTIMATION,
// errorfunction_linear_tolerance_is_decision_in_T_ESTIMATION,
// errorfunction_linear_use_peirce_criterion,
// errorfunction_linear_tolerance_is_biweight_midvariance,
// errorfunction_linear_tolerance_is_percentagebased_midvariance



// The function returns false on error and true if it is successful.




FOCUSINTERPOLATION_API bool focusposition_Regression(vector<long> x, vector<double> y, long* focpos, double* main_error = NULL, double* main_slope = NULL, double* main_intercept = NULL,
	vector<size_t>* indices_of_used_points = NULL,
	vector<double>* usedpoints_line_x = NULL, vector<double>* usedpoints_line_y = NULL, vector<size_t>* indices_of_removedpoints = NULL, vector<double>* removedpoints_line_x = NULL, vector<double>* removedpoints_line_y = NULL,
	double stop_after_seconds = 60, size_t stop_after_numberofiterations_without_improvement = 2000000, long backslash = 0, double scale = 1.5, bool use_median_regression = false,
	size_t maximum_number_of_outliers = 3, outlier_criterion rejection_method = Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION, double tolerance = 3);



// this class contains code that is used for analyzing the images.


class FOCUSINTERPOLATION_API image
{
public:
	// constructs the image class from a fits file given by filename. The focuser position of the image is either supplied as a
	// parameter or read from the fits file. The class can also store hfd and fwhm values if supplied. These values are, however, not used
	// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
	// When the constructor is called, one either has to supply the focuser position of the given fits file as an argument or
	// the constructor attempts to  read the focuser position from the fits file.
	//In order for this to work, the fits files need to have the focuser position recorded under either one of the following 
   // Keywords: FOCUSPOS, FOCUSERPOS, FOCUSERPOSITION, FOCUSPOSITION, FOCUSMOTORPOSITION, FOCUSMOTORPOS
	image(string* filename, long focuser_position = LONG_MIN, double hfd = DBL_MIN, double fwhm = DBL_MIN);

	// constructs the image class from a fits file object. The focuser position of the image is either supplied as a
	// parameter or read from the fits file. The class can also store hfd and fwhm values if supplied. These values are, however, not used
	// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
	// When the constructor is called, one either has to supply the focuser position of the given fits file as an argument or
	// the constructor attempts to  read the focuser position from the fits file.
	//In order for this to work, the fits files need to have the focuser position recorded under either one of the following 
   // Keywords: FOCUSPOS, FOCUSERPOS, FOCUSERPOSITION, FOCUSPOSITION, FOCUSMOTORPOSITION, FOCUSMOTORPOS
	image(fitsfile* fptr, long focuser_position = LONG_MIN, double hfd = DBL_MIN, double fwhm = DBL_MIN);


	//constructs the image class from a double array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
	// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
	image(size_t width, size_t height, long focuser_position, vector <double>* imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN);

	//constructs the image class from a float array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
	// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
	image(size_t width, size_t height, long focuser_position, vector <float>* imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN);

	//constructs the image class from a long long array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
	// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
	image(size_t width, size_t height, long focuser_position, vector <long long>* imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN);

	// constructs the image class from a long array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
	// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
	image(size_t width, size_t height, long focuser_position, vector <long>* imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN);

	// constructs the image class from a short array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
	// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
	image(size_t width, size_t height, long focuser_position, vector <short>* imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN);

	// constructs the image class from an int8_t array. Width and height are the width and height of the image. The focuser position of the image  must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used
	// by the curve fitting procedures. But they may be supplied by applications with an own image analysis.
	image(size_t width, size_t height, long focuser_position, vector <int8_t>* imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN);

	//returns the power of an image in fourier mode.
	double invpower();
	//returns the power of an image in fourier mode.
	double power();
	//returns the focuser position of the image class
	long focuser_position();

	// returns the hfd that can be supplied by the user in the constructor(note that this value is not used in the algorithms, but can be stored as an option, for applications which
	// have their own hdf analysis algorithm.)
	double hfd();

	// returns the fwhm (note that this value is not used in the algorithms, but can be stored as an option, for applications which
	// have their own hdf analysis algorithm.)
	double fwhm();

	//returns the pstatus variable. status=0 means the image class was successfully constructed.
	int status();
private:
	int pstatus;
	long pfocuser_position;
	double phfd;
	double pfwhm;
	double pinvpower;
	double ppower;
	double fouriertransform(vector<float>* p4, vector<double>* p5, size_t dimension1, size_t dimension2,int bitpix);
	double datafilling0(fitsfile* fptr, long focuser_position = LONG_MIN);
};


// focusposition_Regression2 computes the focus point from a vector of validly constructed image classes which must all have status 0.
// these image classes can be constructed from a path to a fits file, a fits file structure, or an array with image data. The constructor
// then computes the power spectrum, with which focusposition_Regression2 can work.
// The parameters of focusposition_Regression2 are similar as in focusposition_Regression.
// 
// However, focusposition_Regression2 fits a function  
// 
// Power=1/ (alpha*(x-focpos)^2+gamma) 
// 
// if the pointer to beta and theta is not zero, this fit is used as an initial value in a Levenberg-Marquardt algorithm
// that fits the function 
// 
// Power=1/ (alpha*(x-focpos)^2+gamma) +exp(-theta/(alpha*(x-focpos)^2+gamma))/(alpha*(x-focpos)^2+gamma)+beta
// 
// if the pointer to theta is zero, the function
// 
// Power=1/ (alpha*(x-focpos)^2+gamma)+beta
// 
// is fitted.
// 
// If the pointer to beta is zero, the function
// 
// Power=1/ (alpha*(x-focpos)^2+gamma) +exp(-theta/(alpha*(x-focpos)^2+gamma))/(alpha*(x-focpos)^2+gamma)
// 
// is fitted.
// 
// Extensive information about the Levenberg - Marquardt algorithm and some techniques used by this application in order to improve convergence can be found in
// Transtrum, Mark K; Sethna, James P(2012). "Improvements to the Levenberg-Marquardt algorithm for nonlinear least-squares minimization".arXiv:1201.5885



FOCUSINTERPOLATION_API bool focusposition_Regression2(std::vector<image>* images, long* focpos, double* main_error=NULL, double* alpha=NULL, double*beta=NULL, double* gamma=NULL, double* theta = NULL,
	vector<size_t>* indices_of_used_points=NULL, vector<size_t>* indices_of_removedpoints=NULL, double stop_after_seconds=60, size_t stop_after_numberofiterations_without_improvement=2000000, long backslash=0, double scale=1, bool use_median_regression=false,
	size_t maximum_number_of_outliers=3, outlier_criterion rejection_method = Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION, double tolerance=3);
// The function findbackslash_Regression finds the focuser backslash from two measurements of the best focus positions. It returns true if successful.
// The idea to correct for the backslash in this way was first published by Jim Hunt (username JST200) at https://aptforum.com/phpbb/viewtopic.php?p=26265#p26265

// The function expects 2 data sets of motor positions x1 and x2 and hfd values y1 and y2. The function then calls focusposition_Regression and makes a curve fit for both sets.
// If the function is successful, the variable value backslash contains the focuser backslash given by the difference 
//
// optimal_focus_point_2 - optimal_focus_point1

// the parameters alpha1, beta1,gamma1, indicesofusedpoints1, used_points1_line_x, used_points1_line_y, indicesofremovedpoints1, removedpoints1_line_x, removedpoints1_line_y
// and alpha2, beta2, gamma2, indicesofusedpoints2, used_points2_line_x, used_points2_line_y, indicesofremovedpoints2, removedpoints2_line_x, removedpoints2_line_y
// have the same meaning as the corresponding return parameters of focusposition_Regression2, just that they are for the two separate datasets x1,y1 and x2,y2.

// The parameters stop_after_seconds, stop_after_numberofiterations_without_improvement, maximum_number_of_outliers, tolerance, scale, use_median_regression, use_median_regression
// have the same meaning as the corresponding parameters in focusposition_Regression and are used for the fits of both datasets.

FOCUSINTERPOLATION_API bool findbackslash_Regression(long* backslash,
	vector<long> x1, vector<double> y1, vector<long> x2, vector<double> y2, double* main_error1=NULL, double* main_slope1=NULL, double* main_intercept1=NULL, vector<size_t>*indicesofusedpoints1=NULL,
	vector<double>*used_points1_line_x=NULL, vector<double>*used_points1_line_y=NULL, vector<size_t>*indicesofremovedpoints1=NULL, vector<double>*removedpoints1_line_x=NULL, vector<double>*removedpoints1_line_y=NULL,
	double* main_error2=NULL, double* main_slope2=NULL, double* main_intercept2=NULL, vector<size_t>*indicesofusedpoints2 = NULL,
	vector<double>*used_points2_line_x=NULL, vector<double>*used_points2_line_y=NULL, vector<size_t>*indicesofremovedpoints2=NULL, vector<double>*removedpoints2_line_x=NULL, vector<double>*removedpoints2_line_y=NULL,
	double stop_after_seconds = 60, size_t stop_after_numberofiterations_without_improvement = 2000000, double scale = 1.5 ,bool use_median_regression = false,
	size_t maximum_number_of_outliers = 3, outlier_criterion rejection_method = Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION, double tolerance = 3);


// The function findbackslash_Regression2 finds the focuser backslash from two interpolations of the best focus positions. In contrast to findbackslash_Regression, it accepts 2 vectors of image classes from 2 consecutive measurements where the 
// focus motor moved in different directions. The parameter theta is currently not used.


FOCUSINTERPOLATION_API bool findbackslash_Regression2(long* backslash,
	vector<image>*images1, vector<image>*images2, double* main_error1 = NULL, double* alpha1 = NULL, double *beta1=NULL,double* gamma1 = NULL, double*theta1=NULL,vector<size_t>*indicesofusedpoints1 = NULL, vector<size_t>*indicesofremovedpoints1 = NULL, 
	double* main_error2 = NULL, double* alpha2 = NULL, double*beta2=NULL, double* gamma2 = NULL, double* theta2=NULL, vector<size_t>*indicesofusedpoints2 = NULL, vector<size_t>*indicesofremovedpoints2 = NULL, 
	double stop_after_seconds = 60, size_t stop_after_numberofiterations_without_improvement = 2000000, double scale = 1.0, bool use_median_regression = false,
	size_t maximum_number_of_outliers = 3, outlier_criterion rejection_method = Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION, double tolerance = 3);
