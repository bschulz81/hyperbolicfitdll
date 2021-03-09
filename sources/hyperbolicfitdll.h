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

#pragma once

#ifdef APTHYPERBOLICFITDLL_EXPORTS
#define HYPERBOLICFIT_API __declspec(dllexport)
#else
#define HYPERBOLICFIT_API __declspec(dllimport)
#endif

using namespace std;
// The function focusposition_Regression interpolates the optimal focuser position from a fit with symmetric hyperbolas based on a RANSAC algorithm that utilizes a linear regression with a least square 
// error comparison. If specified, the RANSAC can also use repeaded median regression.
// The algorithm for the linear regression that is used by the RANSAC was first published by Stephen King (username STEVE333), at https://aptforum.com/phpbb/viewtopic.php?p=25998#p25998

// focusposition_Regression expects the following arguments:

// a vector <long> x with motor positions at which the hfd number was measured. 
// a vector <double> y. The hfd value of a star or a star field corresponding to motor position x [ i ] is assumed to be y [ i ] .
// both x and y should have at least 4 points.

// If the function terminates successfully, the pointer focpos contains the estimated focusposition as a variable of type long that the function computes. This pointer must not be NULL if the
// function is to terminate successfully

// main_error is a pointer to a double value. If the pointer is not NULL, the value of main_error will be the sum of the absolute values of the differences between the best fit and the measurement data. 


// main_slope and main_intercept are pointers to double values. If the pointers are not NULL, their values will be the slope and intercept of a line given by
// Y=slope*(x-focpos)^2+intercept. The square-root f=sqrt(Y) is the fitted hyperbola.

// indices_of_used_points is a pointer to a vector which, if not NULL, will contain the indices i of the points in x [i] and y [i] that were used for the fit.

// usedpoints_line_x and usedpoints_line_y are pointers to two vectors for the datatype double. If the pointers are not NULL, the vectors contain the points which are used in the fit
// in a coordinate system where the hyperbola is a line. One can plot them together with the line  given by Y=slope*(x-focpos)^2+intercept

// indices_of_removedpoints is a pointer to a vector which, if not NULL, will contain the indices i of the points in x [i] and y [i] that were not used for the fit

// removedpoints_line_x and removedpoints_line_y are pointers to two vectors. If the pointers are not NULL, the vectors contain the points which were not used in the fit
// in the coordinate system where the hyperbola is a line. One can plot them together with the line  given by Y=slope*(x-focpos)^2+intercept


// The algorithm works by first selecting combination of points and fitting them to a hyperbola. This initial hyperbola is then corrected with contributions from other points. A point outside a combination is added
// if its error from the initial fit is deemed not to be an outlier based on various statistical methods. A new fitt with the added point is then made and the process is repeated with another initial combination of points.
// 
// The initial combination is selected randomly if the binominal coefficient of the number of points over the number of outliers is larger than 20 over 10. Otherwise, the combinations are searched deterministically.
// 
// stop_after_seconds is a parameter that stops the RANSAC after a given time in seconds has elapsed.
// stop_after_numberofiterations_without_improvement is a parameter that lets the RANSAC stop after it has iterated by stop_after_numberofiterations_without_improvement iterations
// without a further improvement of the error. Note that this parameter is not the iteration number, but it is the number of iterations without further improvement.
// 
// The parameters stop_after_seconds and stop_after_numberofiterations_without_improvement are only used if the binominal coefficient of the number of points over the number of outliers is larger than 20 over 10.


// backslash is a parameter that can contain the focuser backslash in steps. The best focus position is corrected with respect to this backslash. If you already have taken account of
// the focuser backslash, for example by setting a suitable overshoor or a final_inwards_movement in APT or a different software or hardware correction of the backslash, set this parameter to 0


// scale is a parameter of the type double that specifies the size of the interval of motor positions where the best focusposition is searched for. 
// The default of scale is 1.
// Sometimes, the focus point may be outside of the interval of motor positions where the hfd was measured.
// let  middle =(max + min) / 2 and max and min be the maximum and minimum motorposition where a hfd was measured.
// If an initial search finds the best focus point exactly at the right edge of the measurement interval of motor positions,
// then, if scale>1, the right side of the interval where the best focus is searched is enlarged. The new right boundary is then given by
// max = middle + (max - middle) * abs(scale)
// Similarly, if an initial search finds the best focus point exactly on the left side of the measurement interval, the search interval is enlarged, with the new left boundary given by
// min = (middle - (middle - min) * abs(scale).

// use_median_regression is a parameter that specifies whether the RANSAC uses a simple linear regression or a median regression.
// Repeated median regression is slightly more stable against small outliers if one does not use the RANSAC algorithm.

// maximum_number_of_outliers is a parameter that specifies how many outliers the ransac can maximally throw away.


// rejection_method  is a parameter that specifies the method which is used to reject outliers.
// Assume you have n datapoints. The algorithm works by searching through either all or (if the binominal coefficient of points over the number of outliers is larger than 20 over 10) randomly generated so - called minimal combinations
// of m=n - maximum_number_of_outliers points.
// 
// 
// The algorithm searches for the best combination of points with the lowest error, based on linear regression, or repeated median regression.
// For each minimal combination, the points outside of this minimal set of m points are considered. 

// The error between the fit w of a minimal combination and a measurement at a motor position x is given by err_p=p(x)-w(x). 

// If

// rejection_method==no_rejection, then the function uses every point for the fit.


// rejection_method==tolerance_is_maximum_squared_error, 

// then a point p outside of the  minimal combination are only added to the final fit if its squared error err_p*err_p fulfills

// err_p*err_p<=abs(tolerance).

// The RANSAC then computes an overall fit with this combination, and then constructs a new set based on a different minimal mode, until the best combination of points was found. 

// To specify the tolerable error directly is useful if the largest tolerable error is known, e.g. from a measurement of the seeing. By taking a series of images, 
// an application may measure the random deviation of the hfd from the average that stars have. With outlier detection methods, one can remove the outliers of large amplitudes.
// setting tolerance_in_sigma_units=false, one can supply a // maximally tolerable deviation from the average hfd directly to the library as the tolerance parameter.
// The tolerance value then corresponds to the absolute value of the maximally tolerable hfd deviation.

// If 
// rejection_method== tolerance_multiplies_standard_deviation_of_error,

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

// rejection_method==tolerance_is_decision_in_MAD_ESTIMATION, or 
// rejection_method==tolerance_is_decision_in_S_ESTIMATION, or  
// rejection_method==tolerance_is_decision_in_Q_ESTIMATION, or
// rejection_method==tolerance_is_decision_in_T_ESTIMATION,

// then, MAD, S, Q or T estimators are used. 

// A point is then added to a minimal combination if

// abs(err_p-median(err_p))/ estimator <= tolerance

// where median(err_p) is the median of the errors, and estimator is then the MAD, Q, S or T estimator. Q and S estimators are better suited for asymmetric distributions.

// if MAD, Q, or S estimators are used. the tolerance parameter should be around 2...3.


// if 
// rejection_method==tolerance_is_significance_in_Grubbs_test,
//
// then a point p is added to a minimal combination if its error err_p is not an outlier according to the Grubbs test.
// The tolerance parameter is here a significance value based on student's distribution.
// This means, in order to remove an outlier with 90% significance with Grubb's method, one should set tolerance=0.1, for 80% significance, one should select 0.2.

// if 
// rejection_method==use_peirce_criterion, 

// then a point is added to a minimal combination if it does fulfill Peirce's criterion for not being an outlier. 

// For this the rejection method according to Peirce, the tolerance parameter is not used. 
// Note that number of outliers that Peirce's criterion finds is often strongly influenced by the parameter maximum_number_of_outliers

// tolerance is a parameter that is used for the different rejection methods.
 
// The function returns false on error and true if it is successful.

typedef enum
{
	no_rejection,
	tolerance_is_maximum_squared_error,
	tolerance_multiplies_standard_deviation_of_error,
	tolerance_is_significance_in_Grubbs_test,
	tolerance_is_decision_in_MAD_ESTIMATION,
	tolerance_is_decision_in_S_ESTIMATION,
	tolerance_is_decision_in_Q_ESTIMATION,
	tolerance_is_decision_in_T_ESTIMATION,
	use_peirce_criterion
}outlier_criterion;

HYPERBOLICFIT_API bool focusposition_Regression(vector<long> x, vector<double> y, long* focpos, double* main_error = NULL, double* main_slope = NULL, double* main_intercept = NULL,
	vector<size_t>* indices_of_used_points = NULL,
	vector<double>* usedpoints_line_x = NULL, vector<double>* usedpoints_line_y = NULL, vector<size_t>* indices_of_removedpoints = NULL, vector<double>* removedpoints_line_x = NULL, vector<double>* removedpoints_line_y = NULL,
	double stop_after_seconds = 30, size_t stop_after_numberofiterations_without_improvement = 2000, long backslash = 0, double scale = 1.5, bool use_median_regression = false, 
	size_t maximum_number_of_outliers = 3, outlier_criterion rejection_method= tolerance_multiplies_standard_deviation_of_error, double tolerance = 3);

// The function findbackslash finds the focuser backslash from two measurements of the best focus positions. It returns true if successful.
// The idea to correct for the backslash in this way was first published by Jim Hunt (username JST200) at https://aptforum.com/phpbb/viewtopic.php?p=26265#p26265

// The function expects 2 data sets of motor positions x1 and x2 and hfd values y1 and y2. The function then calls focusposition_Regression and makes a curve fit for both sets.
// If the function is successful, the variable value backslash contains the focuser backslash given by the difference 
//
// optimal_focus_point_2 - optimal_focus_point1

// the parameters main_slope1, main_intercept1, indicesofusedpoints1, used_points1_line_x, used_points1_line_y, indicesofremovedpoints1, removedpoints1_line_x, removedpoints1_line_y
// and main_slope2, main_intercept2, indicesofusedpoints2, used_points2_line_x, used_points2_line_y, indicesofremovedpoints2, removedpoints2_line_x, removedpoints2_line_y
// have the same meaning as the corresponding return parameters of focusposition_Regression, just that they are for the two separate datasets x1,y1 and x2,y2.

// The parameters stop_after_seconds, stop_after_numberofiterations_without_improvement, maximum_number_of_outliers, tolerance, scale, use_median_regression, use_median_regression
// have the same meaning as the corresponding parameters in focusposition_Regression and are used for the fits of both datasets.

extern "C" HYPERBOLICFIT_API bool findbackslash_Regression(long* backslash,
	vector<long> x1, vector<double> y1, vector<long> x2, vector<double> y2, double* main_error1=NULL, double* main_slope1=NULL, double* main_intercept1=NULL, vector<size_t>*indicesofusedpoints1=NULL,
	vector<double>*used_points1_line_x=NULL, vector<double>*used_points1_line_y=NULL, vector<size_t>*indicesofremovedpoints1=NULL, vector<double>*removedpoints1_line_x=NULL, vector<double>*removedpoints1_line_y=NULL,
	double* main_error2=NULL, double* main_slope2=NULL, double* main_intercept2=NULL, vector<size_t>*indicesofusedpoints2 = NULL,
	vector<double>*used_points2_line_x=NULL, vector<double>*used_points2_line_y=NULL, vector<size_t>*indicesofremovedpoints2=NULL, vector<double>*removedpoints2_line_x=NULL, vector<double>*removedpoints2_line_y=NULL,
	double stop_after_seconds = 30, size_t stop_after_numberofiterations_without_improvement = 2000, double scale = 1.5, bool use_median_regression = false,
	size_t maximum_number_of_outliers = 3, outlier_criterion rejection_method = tolerance_multiplies_standard_deviation_of_error, double tolerance = 2);