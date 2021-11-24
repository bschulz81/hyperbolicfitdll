// Copyright(c) < 2021 > 
// <Benjamin Schulz> 


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

// This is a test application for the  curve fitting library that can determine the focus point
// of a telescope. The test application expects a folder with fits files at different focuser positions and then interpolates
// the focuser position where the telescope is at focus. 

// In order for this to work, the fits files need to have the focuser position recorded under either one of the following 
// Keywords: FOCUSPOS,FOCUSERPOS,FOCUSERPOSITION,FOCUSPOSITION,FOCUSMOTORPOSITION,FOCUSMOTORPOS

// For example a basic usage of this program may be 
// FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w 
// This would analyze the fits files in the given directory. It will make several curve fits for most of the estimators available in the library and write the results into a file.
// If one does just want to have a single curve fit with a certain estimator, e.g. the S estimator, then one may write
// FM.exe -d Path_To_my_fits_files -info 3 -e lts_S -w
// If one just wants to print out the focus point for the fit with the S estimator and nothing else (which may be useful if this is fed into another rrogram, write
// FM.exe -d Path_To_my_fits_files -info 0 -e lts_S -w
// the above examples use a default maximum number of outliers: It is the pointnumber (or the number of images) -7, this means that there are at least
// 7 points which are assumed not to be outliers. The maximum number of outliers can be specified with the -mo option, e.g
// FM.exe -d Path_To_my_fits_files -info 0 -e lts_S -w -mo 3
// would find at maximum 3 outliers. If this number could reduce the number of usable points below 7, it is set to 0. 
// in general, FM.exe stops with an error if there are not at least 7 images provided.
// Currently FM.exe accepts the following parameters

// -d Directoryname  or --Directory Directoryname, where Directoryname is the path to a folder with fits files
// -w  or --Write_to_file if the application should document its output in a file called Directoryname\\results.txt
// -info number or --Infolevel number, where number ranges from 0 [default] to 3, which sets the level of information that the application 
// prints. 0 means it only prints out the focus point. 3 means it prints the full information
// -m or --Median if the curve fit is done with Siegel's median fitting. The default is false because the median slope is somewhat slow
// -sc number or --Scale number, where scale is a floating point number >=1, that enlarges the area where the focus point is searched. The default is 1
// -bs number or --Backslash number, where number is of type long and should be a previously measured focuser backslash that is subtracted from the
//  estimated focus position. The default is 0. 
// -sec number or --Seconds_for_ransac number  where number is of type double and >0. It specifies the number of seconds after which the ransac is stopped.
//  The default is 60
// -it number or --Iterations_for_ransac number where number is of type long and >0. It specifies the number of iterations without improvement after 
// which the ransac is stopped. The default is 2000000
// -mo number or --Maximum_number_of_outliers number, where number is of type long>=0. It specifies the largest number of outliers that the algorithm 
// can find if it should solve the trimmed least squares problem (the estimators with the lts_ prefix ). The default is the number n of supplied images minus 7. 
// if n<7, then the number of outliers it can find is set to 0. If the binomial coefficient of n and n-8 is larger than 7054320, then the default number 
// of outliers is given by the largest number k where the binomial coefficient of n,k is smaller than 7054320.
// -tol number or --Tolerance number, where number is a value of type float. It sets a cut off value for the estimators. For all estimators the default
//  is 3.0, except for the estumators Significance_in_Grubbs_test, where it is 0.2,Standard_deviation, where it is 2.0 
// -lb indicates whether a Levenberg-Marquardt algorithm should be used to fit the beta parameter for the function 1/(alpha*(x-x0)^2+gamma)+beta
// -lt indicates whether a Levenberg-Marquardt algorithm should be used to fit a second order term i.e the function 1/(alpha*(x-x0)^2+gamma)+exp(-theta/(alpha*(x-x0)^2+gamma))/(alpha*(x-x0)^2+gamma) or 1 / (alpha * (x - x0) ^ 2 + gamma) + exp(-theta / (alpha * (x - x0) ^ 2 + gamma)) / (alpha * (x - x0) ^ 2 + gamma)+beta if -lb is additionally specified
// -e ESTIMATOR or --Estimator ESTIMATOR defines the estimator for the outlier removal procedure. ESTIMATOR can be one of the following variables:
//  NO_REJECTION or 0 This means that no outliers are removed and the tolerance value is ignored. The same behavior can be set by using -mo 0

// ALL_lts This mmeans that the program solves the least trimmed squares problem, where the most optimal point set is searched for. Provided that the maximum number of outliers is >0) several interpolations are made where most of the estimators are triedand the different results noted.For the standard deviation,
// a tolderance of 2 is used, and for the Grubbs test, a tolerance value of 0.2 is used.For the other estimators, the user supplied value or the default 3 is used.If the maximum number of outliers is 0, then only simple linear regressionand the median slope are tested.
// lts_S or 5 This mmeans the least trimmed squares problem is solved where the selection criteria is the S estimator. This is  the default. The tolerance is the cut-off  in |squared_error-median_squared_error|/S(squared_errors)<=tolerance
// lts_MAXIMUM_SQUARED_ERROR or 1 This means the tolerance is the maximum squared error between the fit and the measured points
// lts_STANDARD_DEVIATION or 2 this means that outliers are rejected if their squared_error with respect to the fit is larger than
//  average_squared_error+tolerance*standard_deviation_of_squared_errors
// lts_SIGNIFICANCE_IN_GRUBBS_TEST or 3 This means that the tolerance value defines the significance level in the Grubbs test. 
// lts_MAD or 4. This means that the MAD estimator is used and the tolerance is the cut-off  in |squared_error-average_squared_error|/MAD(squared_errors)<=tolerance
// lts_Q or 6 This means that the Q estimator is used and the tolerance is the cut-off in |squared_error-median_squared_error|/Q(squared_errors)<=tolerance
// lts_T or 7 This means that the T estimator is used and the tolerance is the cut-off in |squared_error-median_squared_error|/T(squared_errors)<=tolerance
// lts_PEIRCE_CRITERION or 8 This means that the Peirce criterion is used
// lts_BIWEIGHT_MIDVARIANCE or 9 This means that the biweight midvariance estimator is used and the tolerance is the cut-off
// in |squared_error-average_squared_error|/biweight-midvariance(squared_errors)<=tolerance;

// ALL_err_vanishing This means several estimators are tested where errors (not points) are removed from the linear regression and Levenberg-Marquardt if they are deemed outliers.
// err_vanishing_MAXIMUM_SQUARED_ERROR or 10 This means errors beyond a squared error given by tolerance are made to vanish
// err_vanishing_STANDARD_DEVIATION or 11 this means errors are set to zero if their squared_error with respect to the fit is larger than average_abs(error)+tolerance*standard_deviation_of_abs(errors)
// err_vanishing_SIGNIFICANCE_IN_GRUBBS_TEST or 12 This means that errors are set to vanish if they fail a Grubs test with a significance level given by tolerance
// err_vanishing_MAD or 13. This means that errors are set to vanish if they are deemed outliers by the MAD estimator, i.e. if they fail to pass |error-median_errors|/MAD(errors)<=tolerance
// err_vanishing_S or 14  This means that errors are set to vanish if they are deemed outliers by the S estimator, i.e. if they fail to pass |error - median_errors|/S(errors) <= tolerance
// err_vanishing_Q or 15 This means that errors are set to vanish if they are deemed outliers by the Q estimator, i.e. if they fail to pass |error-median_errors|/Q(errors)<=tolerance
// err_vanishing_T or 16 This means that errors are set to vanish if they are deemed outliers by the S estimator, i.e. if they fail to pass |error-median_errors|/T(errors)<=tolerance
// err_vanishing_PEIRCE_CRITERION or 17 This means errors vanish if they do not fulfil the Peirce criterion 
// err_vanishing_BIWEIGHT_MIDVARIANCE or 18 This means that errors are set to vanish if they are outliers according to the biweight-midvariance estimator is used and the tolerance is the cut-off
// in |error-median_errors|/biweight-midvariance(errors)<=tolerance

// ALL_err_linear This means several estimators are tested where errors (not points) are only linearly considered as abs(err) in the least squares linear regression
// and Levenberg-Marquardt if they are deemed outliers.
// err_linear_MAXIMUM_SQUARED_ERROR or 19 This means errors beyond a squared error given by tolerance are linearized
// err_linear_STANDARD_DEVIATION or 20 this means errors are are linearized if their squared_error with respect to the fit is larger than average_abs(error)+tolerance*standard_deviation_of_abs(errors)
// err_linear_SIGNIFICANCE_IN_GRUBBS_TEST or 21 This means that errors are linearized if they fail a Grubs test with a significance level given by tolerance
// err_linear_MAD or 22. This means that errors are set to vanish if they are deemed outliers by the MAD estimator, i.e. if they fail to pass |error-median_errors|/MAD(errors)<=tolerance
// err_linear_S or 23  This means that errors are linearized if they are deemed outliers by the S estimator, i.e. if they fail to pass |error - median_errors|/S(errors) <= tolerance
// err_linear_Q or 24 This means that errors are linearized if they are deemed outliers by the Q estimator, i.e. if they fail to pass |error-median_errors|/Q(errors)<=tolerance
// err_linear_T or 25 This means that errors are linearized if they are deemed outliers by the S estimator, i.e. if they fail to pass |error-median_errors|/T(errors)<=tolerance
// err_linear_PEIRCE_CRITERION or 26 This means errors are linearized if they do not fulfil the Peirce criterion " << endl;
// err_linear_BIWEIGHT_MIDVARIANCE or 27 This means that errors linearized set to vanish if they are outliers according to the biweight-midvariance estimator is used and the tolerance is the cut-off in |error-median_errors|/biweight-midvariance(errors)<=tolerance

// ALL means that all estimators for the ALL_lts,ALL_err_vanishing,ALL_err_linear and NO_REJECTION are tested. The curve fit for NO_REJECTION is done twice, with and without median regression for comparison.

// For example a basic usage of this program may be FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w This would analyze the fits files in the given directory. It will make several curve fits for most of the estimators available in the library and write the results into a file.

#include <filesystem>
#include <iostream>
#include <vector>    
#include <string>
#include <stdlib.h>     
#include <chrono>
#include <fstream> 
#include "focusinterpolation.h"

using namespace std;
namespace fs = std::filesystem;


inline size_t binomial(size_t n, size_t k);
void attempt(vector<image>* imv, double scale, double tolerance, double seconds, ofstream* myfile, bool median, int documentationlevel, size_t maxseconds, long maxiterations, long backslash,bool withbe, bool withtheta,bool mapleformat);

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


void attempt(vector<image>* imv, double scale, size_t outliers, double tolerance, outlier_criterion estimator, ofstream* myfile, bool median, int documentationlevel, double seconds, double maxseconds, long maxiterations, long backslash,bool withbe,bool withtheta,bool mapleformat)
{

	double alpha,gamma;

	double beta,theta;

	vector<double> returnline_x, returnline_y, removedpoints_x, removedpoints_y;
	vector<size_t>usedindices, removedindices;

	double mainerror;
	long focpos;


	auto start = std::chrono::steady_clock::now();

	bool b;

	if (withbe && withtheta)
	{
		b = focusposition_Regression2(imv, &focpos, &mainerror, &alpha, &beta, &gamma, &theta, &usedindices, &removedindices, maxseconds, maxiterations, backslash, scale, median, outliers, estimator, tolerance);
	}
	else if (withtheta && !withbe)
	{
		b = focusposition_Regression2(imv, &focpos, &mainerror, &alpha, NULL, &gamma, &theta, &usedindices, &removedindices, maxseconds, maxiterations, backslash, scale, median, outliers, estimator, tolerance);

	}
	else if (withbe && !withtheta)
	{
		b = focusposition_Regression2(imv, &focpos, &mainerror, &alpha, &beta, &gamma, NULL, &usedindices, &removedindices, maxseconds, maxiterations, backslash, scale, median, outliers, estimator, tolerance);

	}
	else if (!withbe && !withtheta)
	{
		b = focusposition_Regression2(imv, &focpos, &mainerror, &alpha, NULL, &gamma, NULL, &usedindices, &removedindices, maxseconds, maxiterations, backslash, scale, median, outliers, estimator, tolerance);

	}

	if (!b) 

	{
		std::cout << "focusposition2 returned  false						" << endl;
		if (myfile != 0)
		{
			*myfile<< "focusposition2 returned  false						" << endl;
		}
		return;
	}
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	seconds += elapsed_seconds.count();

	if (documentationlevel > 0)
	{
		cout << "focusposition: " << focpos << endl<< endl;
		if (myfile != 0)
		{
			*myfile << "focusposition: " << focpos << endl<< endl;
		}
	}
	else
	{
		cout << focpos ;
		if (myfile != 0)
		{
			*myfile << focpos ;
		}
	}
	if (documentationlevel > 0)
	{
		if (myfile != 0)
		{
			if (!withbe && !withtheta)
			{
				if (!mapleformat)
				{
					*myfile << "Power function parameters y = 1 / (alpha * (x - focusposition) ^ 2 + gamma) " << endl;
					*myfile << "alpha:  " << alpha << ", gamma: " << gamma<< ";" << endl << endl;
				}
				else
				{
					*myfile << "f(x, x0, a, b, c):= 1 / (a * (x - x0) ^ 2 + c);" << endl;
					*myfile << "x0:=" << focpos << "; a:=" << alpha << "; c:=" << gamma << ";" << endl;
				}
			}
			else if (withbe && withtheta)
			{
				if (!mapleformat)
				{
					*myfile << "Power function parameters y = 1 / (alpha * (x - focusposition) ^ 2 + gamma)+exp(-theta/(alpha*(x-focusposition)^2+gamma))/(alpha*(x-focusposition)^2+gamma)+beta  " << endl;
					*myfile << "alpha:  " << alpha << ", beta: " << beta << ", gamma: " << gamma << " , theta: " << theta << ";" << endl << endl;
				}
				else
				{
					*myfile << "f(x, x0, a,b c, d):= 1 / (a * (x - x0) ^ 2 + c) + b+exp(-d/(a * (x - x0) ^ 2 + c))/(a * (x - x0) ^ 2 + c);" << endl;
					*myfile << "x0:=" << focpos << "; a:=" << alpha << "; b:= " << beta << "; c:= " << gamma << "; d:= " << theta << ";" << endl;
				}

			}
			else if (withtheta && !withbe)
			{
				if (!mapleformat)
				{
					*myfile << "Power function parameters y = 1 / (alpha * (x - focusposition) ^ 2 + gamma)+exp(-theta/(alpha*(x-focusposition)^2+gamma))/(alpha*(x-focusposition)^2+gamma)  " << endl;
					*myfile << "alpha:  " << alpha << ", gamma: " << gamma << " , theta: " << theta << ";" << endl << endl;
				}
				else
				{
					*myfile << "f(x, x0, a, c, d):= 1 / (a * (x - x0) ^ 2 + c) + exp(-d/(a * (x - x0) ^ 2 + c))/(a * (x - x0) ^ 2 + c);" << endl;
					*myfile << "x0:=" << focpos << "; a:=" << alpha << "; c:=" << gamma << "; d:= " << theta << ";" << endl;
				}
			}
			else if (withbe && !withtheta)
			{
				if (!mapleformat)
				{
					*myfile << "Power function parameters y = 1 / (alpha * (x - focusposition) ^ 2 + gamma)+beta  " << endl;
					*myfile << "alpha:  " << alpha << ", beta. " << beta << ", gamma: " << gamma << ";" << endl << endl;
				}
				else
				{
					*myfile << "f(x, x0, a, b, c):= 1 / (a * (x - x0) ^ 2 + c) + b;" << endl;
					*myfile << "x0:=" << focpos << "; a:=" << alpha << "; b:=" << beta << "; c:=" << gamma << ";" << endl;
				}
			}
		}
		if (!withbe && !withtheta)
		{
			if (!mapleformat)
			{
				cout << "Power function parameters y = 1 / (alpha * (x - focusposition) ^ 2 + gamma) " << endl;
				cout << "alpha:  " << alpha << ", gamma: " << gamma << ";" << endl << endl;
			}
			else
			{
				cout << "f(x, x0, a, b, c):= 1 / (a * (x - x0) ^ 2 + c);" << endl;
				cout << "x0:=" << focpos << "; a:=" << alpha << "; c:=" << gamma << ";" << endl;
			}
		}
		else if (withbe && withtheta)
		{
			if (!mapleformat)
			{
				cout << "Power function parameters y = 1 / (alpha * (x - focusposition) ^ 2 + gamma)+exp(-theta/(alpha*(x-focusposition)^2+gamma))/(alpha*(x-focusposition)^2+gamma)+beta  " << endl;
				cout << "alpha:  " << alpha << ", beta: " << beta << ", gamma: " << gamma << " , theta: " << theta << ";" << endl << endl;
			}
			else
			{
				cout << "f(x, x0, a,b c, d):= 1 / (a * (x - x0) ^ 2 + c) + b+exp(-d/(a * (x - x0) ^ 2 + c))/(a * (x - x0) ^ 2 + c);" << endl;
				cout << "x0:=" << focpos << "; a:=" << alpha << "; b:= " << beta<< "; c:= " << gamma << "; d:= " << theta << ";" << endl;
			}
			
		}
		else if (withtheta && !withbe)
		{
			if (!mapleformat)
			{
				cout << "Power function parameters y = 1 / (alpha * (x - focusposition) ^ 2 + gamma)+exp(-theta/(alpha*(x-focusposition)^2+gamma))/(alpha*(x-focusposition)^2+gamma)  " << endl;
				cout << "alpha:  " << alpha << ", gamma: " << gamma << " , theta " << theta << ";" << endl << endl;
			}
			else
			{
				cout << "f(x, x0, a, c, d):= 1 / (a * (x - x0) ^ 2 + c) + exp(-d/(a * (x - x0) ^ 2 + c))/(a * (x - x0) ^ 2 + c);" << endl;
				cout << "x0:=" << focpos << "; a:=" << alpha << "; c:=" << gamma << "; d:= "<< theta << ";" << endl;
			}
		}
		else if (withbe && !withtheta)
		{
			if (!mapleformat)
			{
				cout << "Power function parameters y = 1 / (alpha * (x - focusposition) ^ 2 + gamma)+beta  " << endl;
				cout << "alpha:  " << alpha << ", b: " << beta << ", gamma: " << gamma << ";" << endl << endl;
			}
			else
			{
				cout << "f(x, x0, a, b, c):= 1 / (a * (x - x0) ^ 2 + c) + b;" << endl;
				cout << "x0:=" << focpos << "; a:=" << alpha << "; b:=" << beta << "; c:=" << gamma << ";" << endl;
			}

		}

	}
	if (documentationlevel >1)
	{
		cout << "estimated error: " << mainerror << endl<< endl;
		if (myfile != 0)
		{
			*myfile << "estimated error: " << mainerror << endl<< endl;
		}
	}
	if (documentationlevel > 2)
	{
		if (outliers > 0)
		{
			cout << endl << "Outliers" << endl;
			cout << "Motorpositions       ";
			if (myfile != 0)
			{
				*myfile << endl << "Outliers" << endl;
				*myfile << "Motorpositions       ";
			}
			for (size_t i = 0; i < removedindices.size(); i++)
			{

				if (myfile != 0)
				{
					*myfile << (*imv)[removedindices[i]].focuser_position();
					if (i < removedindices.size() - 1)
					{
						*myfile << ", ";
					}
				}
				cout << (*imv)[removedindices[i]].focuser_position();
				if (i < removedindices.size() - 1)
				{
					cout << ", ";
				}

			}
			cout << endl;
			cout << "Power function        ";
			if (myfile != 0)
			{
				*myfile << endl;
				*myfile << "Power function        ";
			}
			for (size_t i = 0; i < removedindices.size(); i++)
			{
				if (myfile != 0)
				{
					*myfile << (*imv)[removedindices[i]].power();
					if (i < removedindices.size() - 1)
					{
						*myfile << ", ";
					}

				}
				cout << (*imv)[removedindices[i]].power();
				if (i < removedindices.size() - 1)
				{
					cout << ", ";
				}

			}
		}

	}
}

int main(int argc, char** argv)
{
	string directorypath;
	long outliers1 = LONG_MIN;
	double tolerance = 3.0;
	outlier_criterion estimator = Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION;

	bool median = false;
	bool writetofile = false;
	int documentationlevel = 0;
	long backslash = 0;
	double maxseconds = 60;
	long maxiterations = 2000000;
	double scale = 1.0;
	bool toleranceset = false;
	bool bwiththeta =false;
	bool bwithbe = false;
	bool btestall_lts = false;
	bool btestall_err_vanish = false;
	bool btestall_err_linear = false;
	bool mapleformat = false;
	string entirecommand;
	if (argc == 1)
	{
		cout << "This is the documentation of the test application of the focus interpolation library. The test application expects the following arguments" << endl;
		cout << "-d Directoryname  or --Directory Directoryname, where Directoryname is the path to a folder with fits files" << endl;
		cout << "-w  or --Write_to_file if the application should document its output in a file called Directoryname\\results.txt" << endl;
		cout << "-info number or --Infolevel number, where number ranges from 0 [default] to 3, which sets the level of information that the application prints. 0 means it only prints out the focus point. 3 means it prints the full information" << endl;
		cout << "-m or --Median if the curve fit is done with Siegel's median fitting. The default is false because the median slope is somewhat slow" << endl;
		cout << "-sc number or --Scale number, where scale is a floating point number >=1, that enlarges the area where the focus point is searched. The default is 1" << endl;
		cout << "-bs number or --Backslash number, where number is of type long and should be a previously measured focuser backslash that is subtracted from the estimated focus position. The default is 0. " << endl;
		cout << "-sec number or --Seconds_for_ransac number  where number is of type double and >0. It specifies the number of seconds after which the ransac is stopped. The default is 60" << endl;
		cout << "-it number or --Iterations_for_ransac number where number is of type long and >0. It specifies the number of iterations without improvement after which the ransac is stopped. The default is 2000000" << endl;
		cout << "-mo number or --Maximum_number_of_outliers number, where number is of type long>=0. It specifies the largest number of outliers that the algorithm can find if it should solve the trimmed least squares problem (the estimators with the lts_ prefix ). The default is the number n of supplied images minus 7. if n<7, then the number of outliers it can find is set to 0. If the binomial coefficient of n and n-8 is larger than 7054320, then the default number of outliers is given by the largest number k where the binomial coefficient of n,k is smaller than 7054320." << endl;
		cout << "-tol number or --Tolerance number, where number is a value of type float. It sets a cut off value for the estimators. For all estimators the default is 3.0, except for the estumators Significance_in_Grubbs_test, where it is 0.2,Standard_deviation, where it is 2.0" << endl;
		cout << "-lb indicates whether a Levenberg-Marquardt algorithm should be used to fit the beta parameter for the function 1/(alpha*(x-x0)^2+gamma)+beta" << endl;
		cout << "-lt indicates whether a Levenberg-Marquardt algorithm should be used to fit a second order term i.e the function 1/(alpha*(x-x0)^2+gamma)+exp(-theta/(alpha*(x-x0)^2+gamma))/(alpha*(x-x0)^2+gamma) or 1 / (alpha * (x - x0) ^ 2 + gamma) + exp(-theta / (alpha * (x - x0) ^ 2 + gamma)) / (alpha * (x - x0) ^ 2 + gamma)+beta if -lb is additionally specified" << endl;
		cout << "-e ESTIMATOR or --Estimator ESTIMATOR defines the estimator for the outlier removal procedure. ESTIMATOR can be one of the following variables:" << endl;
		cout << "NO_REJECTION or 0 This means that no outliers are removed and the tolerance value is ignored. The same behavior can be set by using -mo 0" << endl;

		cout << "ALL_lts This mmeans that the program solves the least trimmed squares problem, where the most optimal point set is searched for. Provided that the maximum number of outliers is >0) several interpolations are made where most of the estimators are triedand the different results noted.For the standard deviation, a tolderance of 2 is used, and for the Grubbs test, a tolerance value of 0.2 is used.For the other estimators, the user supplied value or the default 3 is used.If the maximum number of outliers is 0, then only simple linear regressionand the median slope are tested." << endl;
		cout << "lts_S or 5 This mmeans the least trimmed squares problem is solved where the selection criteria is the S estimator. This is  the default. The tolerance is the cut-off  in |squared_error-median_squared_error|/S(squared_errors)<=tolerance" << endl;
		cout << "lts_MAXIMUM_SQUARED_ERROR or 1 This means the tolerance is the maximum squared error between the fit and the measured points" << endl;
		cout << "lts_STANDARD_DEVIATION or 2 this means that outliers are rejected if their squared_error with respect to the fit is larger than average_squared_error+tolerance*standard_deviation_of_squared_errors" << endl;
		cout << "lts_SIGNIFICANCE_IN_GRUBBS_TEST or 3 This means that the tolerance value defines the significance level in the Grubbs test. " << endl;
		cout << "lts_MAD or 4. This means that the MAD estimator is used and the tolerance is the cut-off  in |squared_error-average_squared_error|/MAD(squared_errors)<=tolerance" << endl;
		cout << "lts_Q or 6 This means that the Q estimator is used and the tolerance is the cut-off in |squared_error-median_squared_error|/Q(squared_errors)<=tolerance" << endl;
		cout << "lts_T or 7 This means that the T estimator is used and the tolerance is the cut-off in |squared_error-median_squared_error|/T(squared_errors)<=tolerance" << endl;
		cout << "lts_PEIRCE_CRITERION or 8 This means that the Peirce criterion is used" << endl;
		cout << "lts_BIWEIGHT_MIDVARIANCE or 9 This means that the biweight midvariance estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/biweight-midvariance(squared_errors)<=tolerance" << endl;
		
		cout << "ALL_err_vanishing This means several estimators are tested where errors (not points) are removed from the linear regression and Levenberg-Marquardt if they are deemed outliers." << endl;
		cout << "err_vanishing_MAXIMUM_SQUARED_ERROR or 10 This means errors beyond a squared error given by tolerance are made to vanish" << endl;
		cout << "err_vanishing_STANDARD_DEVIATION or 11 this means errors are set to zero if their squared_error with respect to the fit is larger than average_abs(error)+tolerance*standard_deviation_of_abs(errors)" << endl;
		cout << "err_vanishing_SIGNIFICANCE_IN_GRUBBS_TEST or 12 This means that errors are set to vanish if they fail a Grubs test with a significance level given by tolerance" << endl;
		cout << "err_vanishing_MAD or 13. This means that errors are set to vanish if they are deemed outliers by the MAD estimator, i.e. if they fail to pass |error-median_errors|/MAD(errors)<=tolerance" << endl;
		cout << "err_vanishing_S or 14  This means that errors are set to vanish if they are deemed outliers by the S estimator, i.e. if they fail to pass |error - median_errors|/S(errors) <= tolerance" << endl;
		cout << "err_vanishing_Q or 15 This means that errors are set to vanish if they are deemed outliers by the Q estimator, i.e. if they fail to pass |error-median_errors|/Q(errors)<=tolerance" << endl;
		cout << "err_vanishing_T or 16 This means that errors are set to vanish if they are deemed outliers by the S estimator, i.e. if they fail to pass |error-median_errors|/T(errors)<=tolerance" << endl;
		cout << "err_vanishing_PEIRCE_CRITERION or 17 This means errors vanish if they do not fulfil the Peirce criterion " << endl;
		cout << "err_vanishing_BIWEIGHT_MIDVARIANCE or 18 This means that errors are set to vanish if they are outliers according to the biweight-midvariance estimator is used and the tolerance is the cut-off in |error-median_errors|/biweight-midvariance(errors)<=tolerance" << endl;
		
		cout << "ALL_err_linear This means several estimators are tested where errors (not points) are only linearly considered as abs(err) in the least squares linear regression and Levenberg-Marquardt if they are deemed outliers." << endl;
		cout << "err_linear_MAXIMUM_SQUARED_ERROR or 19 This means errors beyond a squared error given by tolerance are linearized" << endl;
		cout << "err_linear_STANDARD_DEVIATION or 20 this means errors are are linearized if their squared_error with respect to the fit is larger than average_abs(error)+tolerance*standard_deviation_of_abs(errors)" << endl;
		cout << "err_linear_SIGNIFICANCE_IN_GRUBBS_TEST or 21 This means that errors are linearized if they fail a Grubs test with a significance level given by tolerance" << endl;
		cout << "err_linear_MAD or 22. This means that errors are set to vanish if they are deemed outliers by the MAD estimator, i.e. if they fail to pass |error-median_errors|/MAD(errors)<=tolerance" << endl;
		cout << "err_linear_S or 23  This means that errors are linearized if they are deemed outliers by the S estimator, i.e. if they fail to pass |error - median_errors|/S(errors) <= tolerance" << endl;
		cout << "err_linear_Q or 24 This means that errors are linearized if they are deemed outliers by the Q estimator, i.e. if they fail to pass |error-median_errors|/Q(errors)<=tolerance" << endl;
		cout << "err_linear_T or 25 This means that errors are linearized if they are deemed outliers by the S estimator, i.e. if they fail to pass |error-median_errors|/T(errors)<=tolerance" << endl;
		cout << "err_linear_PEIRCE_CRITERION or 26 This means errors are linearized if they do not fulfil the Peirce criterion " << endl;
		cout << "err_linear_BIWEIGHT_MIDVARIANCE or 27 This means that errors linearized set to vanish if they are outliers according to the biweight-midvariance estimator is used and the tolerance is the cut-off in |error-median_errors|/biweight-midvariance(errors)<=tolerance" << endl;
		
		cout << "ALL means that all estimators for the ALL_lts,ALL_err_vanishing,ALL_err_linear and NO_REJECTION are tested. The curve fit for NO_REJECTION is done twice, with and without median regression for comparison." << endl;

		cout << "For example a basic usage of this program may be FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w This would analyze the fits files in the given directory. It will make several curve fits for most of the estimators available in the library and write the results into a file." << endl;
		cout << "press key to end ";
		string input;
		getline(std::cin, input);
		return 0;
	}
	for (size_t i = 1; i < argc; i++)
	{
		entirecommand.append(argv[i]);
		entirecommand.append(" ");
		if (((string)argv[i] == "-h") || ((string)argv[i] == "--Help"))
		{
			cout << "This is the documentation of the test application of the focus interpolation library. The test application expects the following arguments" << endl;
			cout << "-d Directoryname  or --Directory Directoryname, where Directoryname is the path to a folder with fits files" << endl;
			cout << "-w  or --Write_to_file if the application should document its output in a file called Directoryname\\results.txt" << endl;
			cout << "-info number or --Infolevel number, where number ranges from 0 [default] to 3, which sets the level of information that the application prints. 0 means it only prints out the focus point. 3 means it prints the full information" << endl;
			cout << "-m or --Median if the curve fit is done with Siegel's median fitting. The default is false because the median slope is somewhat slow" << endl;
			cout << "-sc number or --Scale number, where scale is a floating point number >=1, that enlarges the area where the focus point is searched. The default is 1" << endl;
			cout << "-bs number or --Backslash number, where number is of type long and should be a previously measured focuser backslash that is subtracted from the estimated focus position. The default is 0. " << endl;
			cout << "-sec number or --Seconds_for_ransac number  where number is of type double and >0. It specifies the number of seconds after which the ransac is stopped. The default is 60" << endl;
			cout << "-it number or --Iterations_for_ransac number where number is of type long and >0. It specifies the number of iterations without improvement after which the ransac is stopped. The default is 2000000" << endl;
			cout << "-mo number or --Maximum_number_of_outliers number, where number is of type long>=0. It specifies the largest number of outliers that the algorithm can find if it should solve the trimmed least squares problem (the estimators with the lts_ prefix ). The default is the number n of supplied images minus 7. if n<7, then the number of outliers it can find is set to 0. If the binomial coefficient of n and n-8 is larger than 7054320, then the default number of outliers is given by the largest number k where the binomial coefficient of n,k is smaller than 7054320." << endl;
			cout << "-tol number or --Tolerance number, where number is a value of type float. It sets a cut off value for the estimators. For all estimators the default is 3.0, except for the estumators Significance_in_Grubbs_test, where it is 0.2,Standard_deviation, where it is 2.0" << endl;
			cout << "-lb indicates whether a Levenberg-Marquardt algorithm should be used to fit the beta parameter for the function 1/(alpha*(x-x0)^2+gamma)+beta" << endl;
			cout << "-lt indicates whether a Levenberg-Marquardt algorithm should be used to fit a second order term i.e the function 1/(alpha*(x-x0)^2+gamma)+exp(-theta/(alpha*(x-x0)^2+gamma))/(alpha*(x-x0)^2+gamma) or 1 / (alpha * (x - x0) ^ 2 + gamma) + exp(-theta / (alpha * (x - x0) ^ 2 + gamma)) / (alpha * (x - x0) ^ 2 + gamma)+beta if -lb is additionally specified" << endl;
			cout << "-e ESTIMATOR or --Estimator ESTIMATOR defines the estimator for the outlier removal procedure. ESTIMATOR can be one of the following variables:" << endl;
			cout << "NO_REJECTION or 0 This means that no outliers are removed and the tolerance value is ignored. The same behavior can be set by using -mo 0" << endl;
			
			cout << "ALL_lts This mmeans that the program solves the least trimmed squares problem, where the most optimal point set is searched for. Provided that the maximum number of outliers is >0) several interpolations are made where most of the estimators are triedand the different results noted.For the standard deviation, a tolderance of 2 is used, and for the Grubbs test, a tolerance value of 0.2 is used.For the other estimators, the user supplied value or the default 3 is used.If the maximum number of outliers is 0, then only simple linear regressionand the median slope are tested." << endl;
			cout << "lts_S or 5 This mmeans the least trimmed squares problem is solved where the selection criteria is the S estimator. This is  the default. The tolerance is the cut-off  in |squared_error-median_squared_error|/S(squared_errors)<=tolerance" << endl;
			cout << "lts_MAXIMUM_SQUARED_ERROR or 1 This means the tolerance is the maximum squared error between the fit and the measured points" << endl;
			cout << "lts_STANDARD_DEVIATION or 2 this means that outliers are rejected if their squared_error with respect to the fit is larger than average_squared_error+tolerance*standard_deviation_of_squared_errors" << endl;
			cout << "lts_SIGNIFICANCE_IN_GRUBBS_TEST or 3 This means that the tolerance value defines the significance level in the Grubbs test. " << endl;
			cout << "lts_MAD or 4. This means that the MAD estimator is used and the tolerance is the cut-off  in |squared_error-average_squared_error|/MAD(squared_errors)<=tolerance" << endl;
			cout << "lts_Q or 6 This means that the Q estimator is used and the tolerance is the cut-off in |squared_error-median_squared_error|/Q(squared_errors)<=tolerance" << endl;
			cout << "lts_T or 7 This means that the T estimator is used and the tolerance is the cut-off in |squared_error-median_squared_error|/T(squared_errors)<=tolerance" << endl;
			cout << "lts_PEIRCE_CRITERION or 8 This means that the Peirce criterion is used" << endl;
			cout << "lts_BIWEIGHT_MIDVARIANCE or 9 This means that the biweight midvariance estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/biweight-midvariance(squared_errors)<=tolerance" << endl;

			cout << "ALL_err_vanishing This means several estimators are tested where errors (not points) are removed from the linear regression and Levenberg-Marquardt if they are deemed outliers." << endl;
			cout << "err_vanishing_MAXIMUM_SQUARED_ERROR or 10 This means errors beyond a squared error given by tolerance are made to vanish" << endl;
			cout << "err_vanishing_STANDARD_DEVIATION or 11 this means errors are set to zero if their squared_error with respect to the fit is larger than average_abs(error)+tolerance*standard_deviation_of_abs(errors)" << endl;
			cout << "err_vanishing_SIGNIFICANCE_IN_GRUBBS_TEST or 12 This means that errors are set to vanish if they fail a Grubs test with a significance level given by tolerance" << endl;
			cout << "err_vanishing_MAD or 13. This means that errors are set to vanish if they are deemed outliers by the MAD estimator, i.e. if they fail to pass |error-median_errors|/MAD(errors)<=tolerance" << endl;
			cout << "err_vanishing_S or 14  This means that errors are set to vanish if they are deemed outliers by the S estimator, i.e. if they fail to pass |error - median_errors|/S(errors) <= tolerance" << endl;
			cout << "err_vanishing_Q or 15 This means that errors are set to vanish if they are deemed outliers by the Q estimator, i.e. if they fail to pass |error-median_errors|/Q(errors)<=tolerance" << endl;
			cout << "err_vanishing_T or 16 This means that errors are set to vanish if they are deemed outliers by the S estimator, i.e. if they fail to pass |error-median_errors|/T(errors)<=tolerance" << endl;
			cout << "err_vanishing_PEIRCE_CRITERION or 17 This means errors vanish if they do not fulfil the Peirce criterion " << endl;
			cout << "err_vanishing_BIWEIGHT_MIDVARIANCE or 18 This means that errors are set to vanish if they are outliers according to the biweight-midvariance estimator is used and the tolerance is the cut-off in |error-median_errors|/biweight-midvariance(errors)<=tolerance" << endl;

			cout << "ALL_err_linear This means several estimators are tested where errors (not points) are only linearly considered as abs(err) in the least squares linear regression and Levenberg-Marquardt if they are deemed outliers." << endl;
			cout << "err_linear_MAXIMUM_SQUARED_ERROR or 19 This means errors beyond a squared error given by tolerance are linearized" << endl;
			cout << "err_linear_STANDARD_DEVIATION or 20 this means errors are are linearized if their squared_error with respect to the fit is larger than average_abs(error)+tolerance*standard_deviation_of_abs(errors)" << endl;
			cout << "err_linear_SIGNIFICANCE_IN_GRUBBS_TEST or 21 This means that errors are linearized if they fail a Grubs test with a significance level given by tolerance" << endl;
			cout << "err_linear_MAD or 22. This means that errors are set to vanish if they are deemed outliers by the MAD estimator, i.e. if they fail to pass |error-median_errors|/MAD(errors)<=tolerance" << endl;
			cout << "err_linear_S or 23  This means that errors are linearized if they are deemed outliers by the S estimator, i.e. if they fail to pass |error - median_errors|/S(errors) <= tolerance" << endl;
			cout << "err_linear_Q or 24 This means that errors are linearized if they are deemed outliers by the Q estimator, i.e. if they fail to pass |error-median_errors|/Q(errors)<=tolerance" << endl;
			cout << "err_linear_T or 25 This means that errors are linearized if they are deemed outliers by the S estimator, i.e. if they fail to pass |error-median_errors|/T(errors)<=tolerance" << endl;
			cout << "err_linear_PEIRCE_CRITERION or 26 This means errors are linearized if they do not fulfil the Peirce criterion " << endl;
			cout << "err_linear_BIWEIGHT_MIDVARIANCE or 27 This means that errors linearized set to vanish if they are outliers according to the biweight-midvariance estimator is used and the tolerance is the cut-off in |error-median_errors|/biweight-midvariance(errors)<=tolerance" << endl;

			cout << "ALL means that all estimators for the ALL_lts,ALL_err_vanishing,ALL_err_linear and NO_REJECTION are tested. The curve fit for NO_REJECTION is done twice, with and without median regression for comparison." << endl;

			cout << "For example a basic usage of this program may be FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w This would analyze the fits files in the given directory. It will make several curve fits for most of the estimators available in the library and write the results into a file." << endl;
			cout << "press key to end ";
			string input;
			getline(std::cin, input);
			return 0;
			break;
		}
		if (((string)argv[i] == "-w") || ((string)argv[i] == "--Write_to_file"))
		{
			writetofile = true;
		}
		if (((string)argv[i] == "-m") || ((string)argv[i] == "--Median"))
		{
			median = true;
		}
		if ((string)argv[i] == "-lb")
			bwithbe = true;
		if ((string)argv[i] == "-lt")
			bwiththeta = true;
		if ((string)argv[i] == "-mf")
			mapleformat = true;
		if (i < (((size_t)argc) - 1))
		{
			if (((string)argv[i] == "-d") || ((string)argv[i] == "--Directory"))
			{
				directorypath = (string)argv[i + 1];
			}
		
			if (((string)argv[i] == "-sc") || ((string)argv[i] == "--Scale"))
			{
				double d = atof(argv[i + 1]);
				if (d > 1) 
				{
					scale = d;
				}
			}
			
			if (((string)argv[i] == "-tol") || ((string)argv[i] == "--Tolerance"))
			{
				 tolerance = atof(argv[i + 1]);
				 toleranceset = true;
			}
			if (((string)argv[i] == "-bs") || ((string)argv[i] == "--Backslash"))
			{
				long l = atol(argv[i + 1]);
				backslash = l;
			}
			if (((string)argv[i] == "-sec") || ((string)argv[i] == "--Seconds_for_ransac"))
			{
				double l = atof(argv[i + 1]);
				if (l > 0)
				{
					maxseconds = l;
				}
			}
			if (((string)argv[i] == "-it") || ((string)argv[i] == "--Iterations_for_ransac"))
			{
				long l = atol(argv[i + 1]);
				if (l > 0)
				{
					maxiterations = l;
				}
			}
			if (((string)argv[i] == "-info") || ((string)argv[i] == "--Infolevel"))
			{
				int l = atol(argv[i + 1]);
				switch (l)
				{
				case 3:
					documentationlevel = 3;
					break;
				case 2:
					documentationlevel = 2;
					break;
				case 1:
					documentationlevel = 1;
					break;
				case 0:
					documentationlevel = 0;
					break;
				}
			}
			if (((string)argv[i] == "-mo") || ((string)argv[i] == "--Maximum_number_of_outliers"))
			{
				long l = atol(argv[i + 1]);
				if (l > 0)
				{
					outliers1 = l;
				}
			}
			if (((string)argv[i] == "-e") || ((string)argv[i] == "--Estimator"))
			{
				if (((string)argv[i + 1] == "NO_REJECTION") || ((string)argv[i + 1] == "0"))
				{
					estimator = no_rejection;
				}
				if (((string)argv[i + 1] == "lts_MAXIMUM_SQUARED_ERROR") || ((string)argv[i + 1] == "1"))
				{
					estimator = Least_trimmed_squares_tolerance_is_maximum_squared_error;
				}
				if (((string)argv[i + 1] == "lts_STANDARD_DEVIATION") || ((string)argv[i + 1] == "2"))
				{
					estimator = Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error;
				}
				if (((string)argv[i + 1] == "lts_SIGNIFICANCE_IN_GRUBBS_TEST") || ((string)argv[i + 1] == "3"))
				{
					estimator = Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test;
				}
				if (((string)argv[i + 1] == "lts_MAD") || ((string)argv[i + 1] == "4"))
				{
					estimator = Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION;
				}
				if (((string)argv[i + 1] == "lts_S") || ((string)argv[i + 1] == "5"))
				{
					estimator = Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION;
				}
				if (((string)argv[i + 1] == "lts_Q") || ((string)argv[i + 1] == "6"))
				{
				estimator = Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION;
				}
				if (((string)argv[i + 1] == "lts_T") || ((string)argv[i + 1] == "7"))
				{
					estimator = Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION;
				}
				if (((string)argv[i + 1] == "lts_PEIRCE_CRITERION") || ((string)argv[i + 1] == "8"))
				{
					estimator = Least_trimmed_squares_use_peirce_criterion;
				}
				if (((string)argv[i + 1] == "lts_BIWEIGHT_MIDVARIANCE") || ((string)argv[i + 1] == "9"))
				{
					estimator = Least_trimmed_squares_tolerance_is_biweight_midvariance;
				}
				
				if (((string)argv[i + 1] == "ALL"))
				{
					btestall_lts = true;
					btestall_err_vanish = true;
					btestall_err_linear = true;
				}

				if (((string)argv[i + 1] == "ALL_lts"))
				{
					btestall_lts = true;
				}
				if (((string)argv[i + 1] == "ALL_err_vanishing"))
				{
					btestall_err_vanish = true;
				}
				if (((string)argv[i + 1] == "ALL_err_linear"))
				{
					btestall_err_linear = true;
				}
				if (((string)argv[i + 1] == "err_vanishing_maximum_squared_error") || ((string)argv[i + 1] == "10"))
				{
					estimator = errorfunction_vanishing_tolerance_is_maximum_squared_error;
				}
				if (((string)argv[i + 1] == "err_vanishing_STANDARD_DEVIATION") || ((string)argv[i + 1] == "11"))
				{
					estimator = errorfunction_vanishing_tolerance_multiplies_standard_deviation_of_error;
				}
				if (((string)argv[i + 1] == "err_vanishing_SIGNIFICANCE_IN_GRUBBS_TEST") || ((string)argv[i + 1] == "12"))
				{
					estimator = errorfunction_vanishing_tolerance_is_significance_in_Grubbs_test;

				}
				if (((string)argv[i + 1] == "err_vanishing_MAD") || ((string)argv[i + 1] == "13"))
				{
					estimator = errorfunction_vanishing_tolerance_is_decision_in_MAD_ESTIMATION;
				}
				if (((string)argv[i + 1] == "err_vanishing_S") || ((string)argv[i + 1] == "14"))
				{
					estimator = errorfunction_vanishing_tolerance_is_decision_in_S_ESTIMATION;
				}
				if (((string)argv[i + 1] == "err_vanishing_Q") || ((string)argv[i + 1] == "15"))
				{
					estimator = errorfunction_vanishing_tolerance_is_decision_in_Q_ESTIMATION;
				}
				if (((string)argv[i + 1] == "err_vanishing_T") || ((string)argv[i + 1] == "16"))
				{
					estimator = errorfunction_vanishing_tolerance_is_decision_in_T_ESTIMATION;
				}
				if (((string)argv[i + 1] == "err_vanishing_PEIRCE_CRITERION") || ((string)argv[i + 1] == "17"))
				{
					estimator = errorfunction_vanishing_use_peirce_criterion;
				}
				if (((string)argv[i + 1] == "err_vanishing_BIWEIGHT_MIDVARIANCE") || ((string)argv[i + 1] == "18"))
				{
					estimator = errorfunction_vanishing_tolerance_is_biweight_midvariance;
				}
					if (((string)argv[i + 1] == "err_linear_maximum_squared_error") || ((string)argv[i + 1] == "19"))
				{
					estimator = errorfunction_linear_tolerance_is_maximum_squared_error;
				}
				if (((string)argv[i + 1] == "err_linear_STANDARD_DEVIATION") || ((string)argv[i + 1] == "20"))
				{
					estimator = errorfunction_linear_tolerance_multiplies_standard_deviation_of_error;
				}
				if (((string)argv[i + 1] == "err_linear_SIGNIFICANCE_IN_GRUBBS_TEST") || ((string)argv[i + 1] == "21"))
				{
					estimator = errorfunction_linear_tolerance_is_significance_in_Grubbs_test;
				}
				if (((string)argv[i + 1] == "err_linear_MAD") || ((string)argv[i + 1] == "22"))
				{
					estimator = errorfunction_linear_tolerance_is_decision_in_MAD_ESTIMATION;
				}
				if (((string)argv[i + 1] == "err_linear_S") || ((string)argv[i + 1] == "23"))
				{
					estimator = errorfunction_linear_tolerance_is_decision_in_S_ESTIMATION;
				}
				if (((string)argv[i + 1] == "err_linear_Q") || ((string)argv[i + 1] == "24"))
				{
					estimator = errorfunction_linear_tolerance_is_decision_in_Q_ESTIMATION;
				}
				if (((string)argv[i + 1] == "err_linear_T") || ((string)argv[i + 1] == "25"))
				{
					estimator = errorfunction_linear_tolerance_is_decision_in_T_ESTIMATION;
				}
				if (((string)argv[i + 1] == "err_linear_PEIRCE_CRITERION") || ((string)argv[i + 1] == "26"))
				{
					estimator = errorfunction_linear_use_peirce_criterion;
				}
				if (((string)argv[i + 1] == "err_linear_BIWEIGHT_MIDVARIANCE") || ((string)argv[i + 1] == "27"))
				{
					estimator = errorfunction_linear_tolerance_is_biweight_midvariance;
				}
			}
		}
	}


	vector<image>imv;
	double seconds = 0;
	auto start = std::chrono::steady_clock::now();
	long pointnumber = 0;
	if (!fs::exists(directorypath))
	{
		cout << "no valid command line options supplied";
		return -1;
	}

	for (const auto& entry : fs::directory_iterator(directorypath))
	{
		string filename_in = entry.path().string();
		if (filename_in.find(".fit") != string::npos)
		{
			image im(&filename_in);
			if (im.status() == 0)
			{
				imv.push_back(im);
				pointnumber++;
			}
		}
	}
	if (pointnumber== 0)
	{
		cout << "no valid image files could be found";
		return -1;
	}
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	seconds += elapsed_seconds.count();

	


	ofstream myfile1;

	if (writetofile)
	{
		fs::path filename_out = directorypath;
		filename_out /="results.txt";
		myfile1.open(filename_out.c_str());
		if (myfile1.is_open() == false)
		{
			cout << "could not open an output file";
			return -1;
		}

	}


	if (pointnumber < 7)
	{
		cout << "not enough valid image files in the directory";
		if (writetofile)
		{
			myfile1 << "not enough valid image filea in the directory";
		}
		return -1;
	}
	if (estimator < 11)
	{
		if (!((estimator == no_rejection) || (outliers1 == 0)))
		{
			if (outliers1 > pointnumber - 7)
			{
				outliers1 = 0;
			}
			if (outliers1 < 0)
			{
				outliers1 = pointnumber - 7;
			}
			if (outliers1 < 0)
			{
				outliers1 = 0;
				estimator = no_rejection;
			}
			if (binomial(pointnumber, outliers1) > 7054320)
			{
				long k = 0;
				while (binomial(pointnumber, k) <= 7054320)
				{
					k++;
				}
				outliers1 = k;
				if (outliers1 == 0)
				{
					estimator = no_rejection;
				}
			}
		}
	}


	if (documentationlevel > 2)
	{
		cout << "command line was " << entirecommand << endl << endl;
		cout << "Points used for the curve fit:" << endl;
		if (writetofile)
		{
			myfile1 << "command line was " << entirecommand << endl << endl;
			myfile1 << "Points used for the curve fit:" << endl;
		}
		if (mapleformat)
		{
			cout << "[";
		}

		if (writetofile)
			if (mapleformat)
			{
				myfile1 << "[";
			}

		for (size_t i = 0; i < imv.size(); i++)
		{

			if (imv[i].status() == 0)
			{
				if(!mapleformat)
					cout << "Focuser position: " << imv[i].focuser_position() << ", Power function: " << imv[i].power() << endl;
				else 
					cout << "[" << imv[i].focuser_position() << "," <<imv[i].power() << "],";
			}
			if (writetofile)
			{
				if (!mapleformat)
					myfile1 << "Focuser position: " << imv[i].focuser_position() << ", Power function: " << imv[i].power() << endl;
				else
				{
					myfile1 << "[" << imv[i].focuser_position() << "," << imv[i].power();
					if (i + 1 < imv.size())
					{
						myfile1 << "],";
					}
					else 
						myfile1 << "]]";
				}
			}
			else
			{
				cout << "Error. Construction of Image class failed for image " << i << endl;
				if (writetofile)
				{
					myfile1 << "Error. Construction of Image class failed for image " << i << endl;
				}
			}
		}
		cout << endl << endl;
		if (writetofile)
			myfile1 << endl << endl;
	}

	if (btestall_lts || btestall_err_vanish || btestall_err_linear)
	{
		if (writetofile)
		{
			cout << "simple linear regression" << endl;
			if (writetofile) myfile1 << endl << endl << "simple linear regression" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			cout << endl << endl << "Siegel's median slope" << endl;
			myfile1 << endl << endl << "Siegel's median slope" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, &myfile1, true, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta, mapleformat);
		}
		else
		{
			cout << "simple linear regression" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta, mapleformat);


			cout << endl << endl << "Siegel's median slope" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, NULL, true, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta, mapleformat);
		}
	}
	if (btestall_err_vanish)
	{
		if (writetofile)
		{
				cout << endl << endl << "errors removed by MAD estimation" << endl;
				myfile1 << endl << endl << "errors removed by MAD estimation" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_tolerance_is_decision_in_MAD_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);



				cout << endl << endl << "errors removed by S estimator" << endl;
				myfile1 << endl << endl << "errors removed by S estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_tolerance_is_decision_in_S_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "errors removed by Q estimator" << endl;
				myfile1 << endl << endl << "errors removed by Q estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_tolerance_is_decision_in_Q_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);



				cout << endl << endl << "errors removed by T estimator" << endl;
				myfile1 << endl << endl << "errors removed by T estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_tolerance_is_decision_in_T_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "errors removed by biweight Midvariance estimator" << endl;
				myfile1 << endl << endl << "errors removed by biweight Midvariance estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_tolerance_is_biweight_midvariance, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "errors removed by Peirce criterion" << endl;
				myfile1 << endl << endl << "errors removed by Peirce criterion" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_use_peirce_criterion, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "removes errors whose squared error is above average plus 2 standard deviations" << endl;
				myfile1 << endl << endl << "removes errors whose squared error is above average plus 2 standard deviations" << endl;
				attempt(&imv, scale, outliers1, 2, errorfunction_vanishing_tolerance_multiplies_standard_deviation_of_error, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "errors removed by Grubs test at 80% significance" << endl;
				myfile1 << endl << endl << "errors removed by Grubs test at 80% significance" << endl;
				attempt(&imv, scale, outliers1, 0.2, errorfunction_vanishing_tolerance_is_significance_in_Grubbs_test, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

		}
		else
		{
					
				cout << endl << endl << "errors removed by MAD estimation" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_tolerance_is_decision_in_MAD_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "errors removed by S estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_tolerance_is_decision_in_S_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "errors removed by Q estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_tolerance_is_decision_in_Q_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "errors removed by T estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_tolerance_is_decision_in_T_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "errors removed by biweight Midvariance estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_tolerance_is_biweight_midvariance, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				

				cout << endl << endl << "errors removed by Peirce criterion" << endl;
				attempt(&imv, scale, outliers1, tolerance, errorfunction_vanishing_use_peirce_criterion, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "removes errors whose squared error is above average plus 2 standard deviations" << endl;
				attempt(&imv, scale, outliers1, 2, errorfunction_vanishing_tolerance_multiplies_standard_deviation_of_error, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


				cout << endl << endl << "errors removed by Grubs test at 80% significance" << endl;
				attempt(&imv, scale, outliers1, 0.2, errorfunction_vanishing_tolerance_is_significance_in_Grubbs_test, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

		}
	}
	else if (btestall_err_linear)
	{
		if (writetofile)
		{
			cout << endl << endl << "errors removed by MAD estimation" << endl;
			myfile1 << endl << endl << "errors removed by MAD estimation" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_decision_in_MAD_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);



			cout << endl << endl << "errors removed by S estimator" << endl;
			myfile1 << endl << endl << "errors removed by S estimator" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_decision_in_S_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "errors removed by Q estimator" << endl;
			myfile1 << endl << endl << "errors removed by Q estimator" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_decision_in_Q_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);



			cout << endl << endl << "errors removed by T estimator" << endl;
			myfile1 << endl << endl << "errors removed by T estimator" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_decision_in_T_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "errors removed by biweight Midvariance estimator" << endl;
			myfile1 << endl << endl << "errors removed by biweight Midvariance estimator" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_biweight_midvariance, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);




			cout << endl << endl << "errors removed by Peirce criterion" << endl;
			myfile1 << endl << endl << "errors removed by Peirce criterion" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_use_peirce_criterion, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "removes errors whose squared error is above average plus 2 standard deviations" << endl;
			myfile1 << endl << endl << "removes errors whose squared error is above average plus 2 standard deviations" << endl;
			attempt(&imv, scale, outliers1, 2, errorfunction_linear_tolerance_multiplies_standard_deviation_of_error, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "errors removed by Grubs test at 80% significance" << endl;
			myfile1 << endl << endl << "errors removed by Grubs test at 80% significance" << endl;
			attempt(&imv, scale, outliers1, 0.2, errorfunction_linear_tolerance_is_significance_in_Grubbs_test, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

		}
		else
		{

			cout << endl << endl << "errors linearized by MAD estimation" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_decision_in_MAD_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "errors linearized by S estimator" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_decision_in_S_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "errors linearized by Q estimator" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_decision_in_Q_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "errors linearized by T estimator" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_decision_in_T_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "errors linearized by biweight Midvariance estimator" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_biweight_midvariance, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "errors linearized by percentage based midvariance estimator" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_tolerance_is_percentagebased_midvariance, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "errors linearized by Peirce criterion" << endl;
			attempt(&imv, scale, outliers1, tolerance, errorfunction_linear_use_peirce_criterion, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "linearizes errors whose squared error is above average plus 2 standard deviations" << endl;
			attempt(&imv, scale, outliers1, 2, errorfunction_linear_tolerance_multiplies_standard_deviation_of_error, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);


			cout << endl << endl << "errors linearized by Grubs test at 80% significance" << endl;
			attempt(&imv, scale, outliers1, 0.2, errorfunction_linear_tolerance_is_significance_in_Grubbs_test, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

		}
	}
	else if (btestall_lts)
	{
		if (writetofile)
		{
			cout <<"simple linear regression" << endl;
			if (writetofile) myfile1 << endl<< endl <<"simple linear regression" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash,bwithbe,bwiththeta,mapleformat);


		
			cout << endl << endl<< "Siegel's median slope" << endl;
			myfile1 << endl << endl<< "Siegel's median slope" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, &myfile1, true, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			if (outliers1 > 0)
			{
				cout << endl << endl << "MAD estimation" << endl;
				myfile1 << endl << endl << "MAD estimation" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

		

				cout << endl << endl << "S estimator" << endl;
				myfile1 << endl << endl << "S estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

	

				cout << endl << endl << "Q estimator" << endl;
				myfile1 << endl << endl << "Q estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

		

				cout << endl << endl << "T estimator" << endl;
				myfile1 << endl << endl << "T estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

		
				cout << endl << endl << "Biweight Midvariance estimator" << endl;
				myfile1 << endl << endl << "Biweight Midvariance estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_tolerance_is_biweight_midvariance, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			
				
			
				cout << endl << endl << "Peirce criterion" << endl;
				myfile1 << endl << endl << "Peirce criterion" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_use_peirce_criterion, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);
				
			
				cout << endl << endl << "removes points whose squared error is above average plus 2 standard deviations" << endl;
				myfile1 << endl << endl << "removes points whose squared error is above average plus 2 standard deviations" << endl;
				attempt(&imv, scale, outliers1, 2, Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			
				cout << endl << endl << "Grubs test at 80% significance" << endl;
				myfile1 << endl << endl << "Grubs test at 80% significance" << endl;
				attempt(&imv, scale, outliers1, 0.2, Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			}
		}
		else
		{
			cout << "simple linear regression" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

		

			cout << endl << endl << "Siegel's median slope" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, NULL, true, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			if (outliers1 > 0)
			{
				cout << endl << endl << "MAD estimation" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

		
				cout << endl << endl << "S estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			
				cout << endl << endl << "Q estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			
				cout << endl << endl << "T estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			
				cout << endl << endl << "Biweight Midvariance estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_tolerance_is_biweight_midvariance, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			
			
			
				cout << endl << endl << "Peirce criterion" << endl;
				attempt(&imv, scale, outliers1, tolerance, Least_trimmed_squares_use_peirce_criterion, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			
				cout << endl << endl << "removes points whose squared error is above average plus 2 standard deviations" << endl;
				attempt(&imv, scale, outliers1, 2, Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			
				cout << endl << endl << "Grubs test at 80% significance" << endl;
				attempt(&imv, scale, outliers1, 0.2, Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);

			}
		}
	}
	else
	{
		if (estimator == Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test)
		{
			if (toleranceset == false)
				tolerance = 0.9;
		}
		if (estimator == Least_trimmed_squares_tolerance_multiplies_standard_deviation_of_error)
		{
			if (toleranceset == false)
				tolerance = 2;
		}
		if (estimator == errorfunction_vanishing_tolerance_is_significance_in_Grubbs_test)
		{
			if (toleranceset == false)
				tolerance = 0.9;
		}
		if (estimator == errorfunction_vanishing_tolerance_multiplies_standard_deviation_of_error)
		{
			if (toleranceset == false)
				tolerance = 2;
		}
		if (estimator == errorfunction_linear_tolerance_is_significance_in_Grubbs_test)
		{
			if (toleranceset == false)
				tolerance = 0.9;
		}
		if (estimator == errorfunction_linear_tolerance_multiplies_standard_deviation_of_error)
		{
			if (toleranceset == false)
				tolerance = 2;
		}
		if (writetofile)
		{
			attempt(&imv, scale, outliers1, tolerance, estimator, &myfile1, median, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);
		}
		else
		{
			attempt(&imv, scale, outliers1, tolerance, estimator, NULL, median, documentationlevel, seconds, maxseconds, maxiterations, backslash, bwithbe, bwiththeta,mapleformat);
		}
	}
	cout << endl;
	if (writetofile)
	{
		myfile1 << endl;
	}

	if (documentationlevel > 0)
	{
		cout  << "Seconds for the computation: " << seconds << endl;
		if (writetofile)
		{
			myfile1  << "Seconds for the computation: " << seconds << endl;;
		}
	}
	if (writetofile)
	{
		myfile1.close();
	}
	cout << "press key to end ";
	string input;
	getline(std::cin, input);
	

	return 0;
}
