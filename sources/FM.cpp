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
// FM.exe -d Path_To_my_fits_files -info 3 -e S -w
// If one just wants to print out the focus point for the fit with the S estimator and nothing else (which may be useful if this is fed into another rrogram, write
// FM.exe -d Path_To_my_fits_files -info 0 -e S -w
// the above examples use a default maximum number of outliers: It is the pointnumber (or the number of images) -7, this means that there are at least
// 7 points which are assumed not to be outliers. The maximum number of outliers can be specified with the -mo option, e.g
// FM.exe -d Path_To_my_fits_files -info 0 -e S -w -mo 3
// would find at maximum 3 outliers. If this number could reduce the number of usable points below 7, it is set to 0. 
// in general, FM.exe stops with an error if there are not at least 7 images provided.
// Currently FM.exe accepts the following parameters

// -d Directoryname  or --Directory Directoryname, where Directoryname is the path to a folder with fits files
// -w  or --Write_to_file if the application should document its output in a file called Directoryname\\results.txt
// -info number or --Infolevel number, where number ranges from 0 [default] to 3, which sets the level of information that the application prints. 0 means it only prints out the focus point. 3 means it prints the full information
// -m or --Median if the curve fit is done with Siegel's median fitting. The default is false because the median slope is somewhat slow
// -sc number or --Scale number, where scale is a floating point number >=1, that enlarges the area where the focus point is searched. The default is 1
// -bs number or --Backslash number, where number is of type long and should be a previously measured focuser backslash that is subtracted from the estimated focus position. The default is 0. 
// -sec number or --Seconds_for_ransac number  where number is of type double and >0. It specifies the number of seconds after which the ransac is stopped. The default is 60
// -it number or --Iterations_for_ransac number where number is of type long and >0. It specifies the number of iterations without improvement after which the ransac is stopped. The default is 2000000
// -mo number or --Maximum_number_of_outliers number, where number is of type long>=0. It specifies the largest number of outliers that the algorithm can find. The default is the number n of supplied images minus 7. if n<7, then the number of outliers it can find is set to 0. If the binomial coefficient of n and n-8 is larger than 7054320, then the default number of outliers is given by the largest number k where the binomial coefficient of n,k is smaller than 7054320.
// -tol number or --Tolerance number, where number is a value of type float. It sets a cut off value for the estimators. For all estimators the default is 3.0, except for the estumators Significance_in_Grubbs_test, where it is 0.2,Standard_deviation, where it is 2.0
// -e ESTIMATOR or --Estimator ESTIMATOR defines the estimator for the outlier removal procedure. ESTIMATOR can be one of the following variables:
// S or 5 This mmeans the S estimator is used, which is the default. The tolerance is the cut-off  in |squared_error-average_squared_error|/S(squared_errors)<=tolerance
// ALL This mmeans that, provided the maximum number of outliers is >0) several interpolations are made where most of the estimators are tried and the different results noted. For the standard deviation, a tolderance of 2 is used, and for the Grubbs test, a tolerance value of 0.2 is used. For the other estimators, the user supplied value or the default 3 is used. If the maximum number of outliers is 0, then only simple linear regression and the median slope are tested.
// NO_REJECTION or 0 This means that no outliers are removed and the tolerance value is ignored. The same behavior can be set by using -mo 0
// MAXIMUM_SQUARED_ERROR or 1 This means the tolerance is the maximum squared error between the fit and the measured points
// STANDARD_DEVIATION or 2 this means that outliers are rejected if their squared_error with respect to the fit is larger than average_squared_error+tolerance*standard_deviation_of_squared_errors
// SIGNIFICANCE_IN_GRUBBS_TEST or 3 This means that the tolerance value defines the significance level in the Grubbs test. 
// MAD or 4. This means that the MAD estimator is used and the tolerance is the cut-off  in |squared_error-average_squared_error|/MAD(squared_errors)<=tolerance
// Q or 6 This means that the Q estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/Q(squared_errors)<=tolerance
// T or 7 This means that the T estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/T(squared_errors)<=tolerance
// PEIRCE_CRITERION or 8 This means that the Peirce criterion is used
// BIWEIGHT_MIDVARIANCE or 9 This means that the biweight midvariance estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/biweight-midvariance(squared_errors)<=tolerance
// PERCENTAGE_BASED_MIDVARIANCE or 10 This means that the percentage-based midvariance estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/percentage_based_midvariance(squared_errors)<=tolerance
// For example a basic usage of this program may be FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w This would analyze the fits files in the given directory. It will make several curve fits for most of the estimators available in the library and write the results into a file.


#include<filesystem>
#include <iostream>
#include <vector>    
#include <string>
#include <stdlib.h>     
#include <chrono>
#include <fstream> 
#include "hyperbolicfitdll.h"

using namespace std;
namespace fs = std::experimental::filesystem;


inline size_t binominal(size_t n, size_t k);
void attempt(vector<image>* imv, double scale, double tolerance, double seconds, ofstream* myfile, bool median, int documentationlevel, size_t maxseconds, long maxiterations, long backslash);

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




void attempt(vector<image>* imv, double scale, size_t outliers, double tolerance, outlier_criterion estimator, ofstream* myfile, bool median, int documentationlevel, double seconds, double maxseconds, long maxiterations, long backslash) 
{

	double mainslope, mainintercept;
	vector<double> returnline_x, returnline_y, removedpoints_x, removedpoints_y;
	vector<size_t>usedindices, removedindices;

	double mainerror;
	long focpos;


	auto start = std::chrono::steady_clock::now();


	if (focusposition_Regression2(imv, &focpos, &mainerror, &mainslope, &mainintercept, NULL, &usedindices, &removedindices, maxseconds, maxiterations, backslash, scale, median, outliers, estimator, tolerance) == false)
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
			*myfile << "Parabola parameters y=slope (x-bestfocus)^2+intercept" << endl;
			*myfile << "intercept:  " << mainintercept << ", slope: " << mainslope << endl << endl;
		}
		cout << "Parabola parameters y=slope (x-bestfocus)^2+intercept" << endl;
		cout << "intercept:  " << mainintercept << ", slope: " << mainslope << endl << endl;
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
		cout << "Points used for the curve fit:"  << endl;
		if (myfile != 0)
		{
			*myfile << "Points used for the curve fit:"  << endl;
		}
		for (size_t i = 0; i < imv->size(); i++)
		{
			if ((*imv)[i].status() == 0)
			{
				cout << "Focuser position: " << (*imv)[i].focuser_position() << ", inverse Power function: " << (*imv)[i].invpower() << endl;
			}
			if (myfile != 0)
			{
				*myfile << "Focuser position: " << (*imv)[i].focuser_position() << ", inverse Power function: " << (*imv)[i].invpower() << endl;
			}
			else
			{
				cout << "Error. Construction of Image class failed for image " << i << endl;
				if (myfile != 0)
				{
					*myfile << "Error. Construction of Image class failed for image " << i << endl;
				}
			}
		}
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
			cout << "inverse Power function        ";
			if (myfile != 0)
			{
				*myfile << endl;
				*myfile << "inverse Power function        ";
			}
			for (size_t i = 0; i < removedindices.size(); i++)
			{
				if (myfile != 0)
				{
					*myfile << (*imv)[removedindices[i]].invpower();
					if (i < removedindices.size() - 1)
					{
						*myfile << ", ";
					}

				}
				cout << (*imv)[removedindices[i]].invpower();
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
	long outliers1 = -LONG_MIN;
	double tolerance = 3.0;
	outlier_criterion estimator = tolerance_is_decision_in_S_ESTIMATION;

	bool testallestimators = false;
	bool median = false;
	bool writetofile = false;
	int documentationlevel = 0;
	long backslash = 0;
	double maxseconds = 60;
	long maxiterations = 2000000;
	double scale = 1.0;
	bool toleranceset = false;
	if (argc == 1)
	{
		cout << "not enough command line arguments"<< endl;
		cout << "This is the documentation of the test application of the focus interpolation library. The test application expects the following arguments" << endl;
		cout << "-d Directoryname  or --Directory Directoryname, where Directoryname is the path to a folder with fits files" << endl;
		cout << "-w  or --Write_to_file if the application should document its output in a file called Directoryname\\results.txt" << endl;
		cout << "-info number or --Infolevel number, where number ranges from 0 [default] to 3, which sets the level of information that the application prints. 0 means it only prints out the focus point. 3 means it prints the full information" << endl;
		cout << "-m or --Median if the curve fit is done with Siegel's median fitting. The default is false because the median slope is somewhat slow" << endl;
		cout << "-sc number or --Scale number, where scale is a floating point number >=1, that enlarges the area where the focus point is searched. The default is 1" << endl;
		cout << "-bs number or --Backslash number, where number is of type long and should be a previously measured focuser backslash that is subtracted from the estimated focus position. The default is 0. " << endl;
		cout << "-sec number or --Seconds_for_ransac number  where number is of type double and >0. It specifies the number of seconds after which the ransac is stopped. The default is 60" << endl;
		cout << "-it number or --Iterations_for_ransac number where number is of type long and >0. It specifies the number of iterations without improvement after which the ransac is stopped. The default is 2000000" << endl;
		cout << "-mo number or --Maximum_number_of_outliers number, where number is of type long>=0. It specifies the largest number of outliers that the algorithm can find. The default is the number n of supplied images minus 7. if n<7, then the number of outliers it can find is set to 0. If the binomial coefficient of n and n-8 is larger than 7054320, then the default number of outliers is given by the largest number k where the binomial coefficient of n,k is smaller than 7054320." << endl;
		cout << "-tol number or --Tolerance number, where number is a value of type float. It sets a cut off value for the estimators. For all estimators the default is 3.0, except for the estumators Significance_in_Grubbs_test, where it is 0.2,Standard_deviation, where it is 2.0" << endl;
		cout << "-e ESTIMATOR or --Estimator ESTIMATOR defines the estimator for the outlier removal procedure. ESTIMATOR can be one of the following variables:" << endl;
		cout << "S or 5 This mmeans the S estimator is used, which is the default. The tolerance is the cut-off  in |squared_error-average_squared_error|/S(squared_errors)<=tolerance" << endl;
		cout << "ALL This mmeans that, provided the maximum number of outliers is >0) several interpolations are made where most of the estimators are tried and the different results noted. For the standard deviation, a tolderance of 2 is used, and for the Grubbs test, a tolerance value of 0.2 is used. For the other estimators, the user supplied value or the default 3 is used. If the maximum number of outliers is 0, then only simple linear regression and the median slope are tested." << endl;
		cout << "NO_REJECTION or 0 This means that no outliers are removed and the tolerance value is ignored. The same behavior can be set by using -mo 0" << endl;
		cout << "MAXIMUM_SQUARED_ERROR or 1 This means the tolerance is the maximum squared error between the fit and the measured points" << endl;
		cout << "STANDARD_DEVIATION or 2 this means that outliers are rejected if their squared_error with respect to the fit is larger than average_squared_error+tolerance*standard_deviation_of_squared_errors" << endl;
		cout << "SIGNIFICANCE_IN_GRUBBS_TEST or 3 This means that the tolerance value defines the significance level in the Grubbs test. " << endl;
		cout << "MAD or 4. This means that the MAD estimator is used and the tolerance is the cut-off  in |squared_error-average_squared_error|/MAD(squared_errors)<=tolerance" << endl;
		cout << "Q or 6 This means that the Q estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/Q(squared_errors)<=tolerance" << endl;
		cout << "T or 7 This means that the T estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/T(squared_errors)<=tolerance" << endl;
		cout << "PEIRCE_CRITERION or 8 This means that the Peirce criterion is used" << endl;
		cout << "BIWEIGHT_MIDVARIANCE or 9 This means that the biweight midvariance estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/biweight-midvariance(squared_errors)<=tolerance" << endl;
		cout << "PERCENTAGE_BASED_MIDVARIANCE or 10 This means that the percentage-based midvariance estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/percentage_based_midvariance(squared_errors)<=tolerance" << endl;
		cout << "For example a basic usage of this program may be FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w This would analyze the fits files in the given directory. It will make several curve fits for most of the estimators available in the library and write the results into a file." << endl;
		return 0;
	}
	for (size_t i = 1; i < argc; i++)
	{
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
			cout << "-mo number or --Maximum_number_of_outliers number, where number is of type long>=0. It specifies the largest number of outliers that the algorithm can find. The default is the number n of supplied images minus 7. if n<7, then the number of outliers it can find is set to 0. If the binomial coefficient of n and n-8 is larger than 7054320, then the default number of outliers is given by the largest number k where the binomial coefficient of n,k is smaller than 7054320." << endl;
			cout << "-tol number or --Tolerance number, where number is a value of type float. It sets a cut off value for the estimators. For all estimators the default is 3.0, except for the estumators Significance_in_Grubbs_test, where it is 0.2,Standard_deviation, where it is 2.0" << endl;
			cout << "-e ESTIMATOR or --Estimator ESTIMATOR defines the estimator for the outlier removal procedure. ESTIMATOR can be one of the following variables:" << endl;
			cout << "S or 5 This mmeans the S estimator is used, which is the default. The tolerance is the cut-off  in |squared_error-average_squared_error|/S(squared_errors)<=tolerance" << endl;
			cout << "ALL This mmeans that, provided the maximum number of outliers is >0) several interpolations are made where most of the estimators are tried and the different results noted. For the standard deviation, a tolderance of 2 is used, and for the Grubbs test, a tolerance value of 0.2 is used. For the other estimators, the user supplied value or the default 3 is used. If the maximum number of outliers is 0, then only simple linear regression and the median slope are tested." << endl;
			cout << "NO_REJECTION or 0 This means that no outliers are removed and the tolerance value is ignored. The same behavior can be set by using -mo 0" << endl;
			cout << "MAXIMUM_SQUARED_ERROR or 1 This means the tolerance is the maximum squared error between the fit and the measured points" << endl;
			cout << "STANDARD_DEVIATION or 2 this means that outliers are rejected if their squared_error with respect to the fit is larger than average_squared_error+tolerance*standard_deviation_of_squared_errors" << endl;
			cout << "SIGNIFICANCE_IN_GRUBBS_TEST or 3 This means that the tolerance value defines the significance level in the Grubbs test. " << endl;
			cout << "MAD or 4. This means that the MAD estimator is used and the tolerance is the cut-off  in |squared_error-average_squared_error|/MAD(squared_errors)<=tolerance" << endl;
			cout << "Q or 6 This means that the Q estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/Q(squared_errors)<=tolerance" << endl;
			cout << "T or 7 This means that the T estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/T(squared_errors)<=tolerance" << endl;
			cout << "PEIRCE_CRITERION or 8 This means that the Peirce criterion is used" << endl;
			cout << "BIWEIGHT_MIDVARIANCE or 9 This means that the biweight midvariance estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/biweight-midvariance(squared_errors)<=tolerance" << endl;
			cout << "PERCENTAGE_BASED_MIDVARIANCE or 10 This means that the percentage-based midvariance estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/percentage_based_midvariance(squared_errors)<=tolerance" << endl;
			cout << "For example a basic usage of this program may be FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w This would analyze the fits files in the given directory. It will make several curve fits for most of the estimators available in the library and write the results into a file." << endl;
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
				if (((string)argv[i + 1] == "MAXIMUM_SQUARED_ERROR") || ((string)argv[i + 1] == "1"))
				{
					estimator = tolerance_is_maximum_squared_error;
				}
				if (((string)argv[i + 1] == "STANDARD_DEVIATION") || ((string)argv[i + 1] == "2"))
				{
					estimator = tolerance_multiplies_standard_deviation_of_error;
				}
				if (((string)argv[i + 1] == "SIGNIFICANCE_IN_GRUBBS_TEST") || ((string)argv[i + 1] == "3"))
				{
					estimator = tolerance_is_significance_in_Grubbs_test;
				}
				if (((string)argv[i + 1] == "MAD") || ((string)argv[i + 1] == "4"))
				{
					estimator = tolerance_is_decision_in_MAD_ESTIMATION;
				}
				if (((string)argv[i + 1] == "S") || ((string)argv[i + 1] == "5"))
				{
					estimator = tolerance_is_decision_in_S_ESTIMATION;
				}
				if (((string)argv[i + 1] == "Q") || ((string)argv[i + 1] == "6"))
				{
				estimator = tolerance_is_decision_in_Q_ESTIMATION;
				}
				if (((string)argv[i + 1] == "T") || ((string)argv[i + 1] == "7"))
				{
					estimator = tolerance_is_decision_in_T_ESTIMATION;
				}
				if (((string)argv[i + 1] == "PEIRCE_CRITERION") || ((string)argv[i + 1] == "8"))
				{
					estimator = use_peirce_criterion;
				}
				if (((string)argv[i + 1] == "BIWEIGHT_MIDVARIANCE") || ((string)argv[i + 1] == "9"))
				{
					estimator = tolerance_is_biweight_midvariance;
				}
				if (((string)argv[i + 1] == "PERCENTAGE_BASED_MIDVARIANCE") || ((string)argv[i + 1] == "10"))
				{
					estimator = tolerance_is_percentagebased_midvariance;
				}
				if (((string)argv[i + 1] == "ALL"))
				{
					testallestimators = true;
				}
			}
		}
	}


	vector<image>imv;
	double seconds = 0;
	auto start = std::chrono::steady_clock::now();
	for (const auto& entry : fs::directory_iterator(directorypath))
	{
		string filename_in = entry.path().string();
		if (filename_in.find(".fit") != string::npos)
		{
			image im(&filename_in);
			imv.push_back(im);
		}
	}

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	seconds += elapsed_seconds.count();

	long pointnumber = 0;
	for (size_t i=0;i< imv.size(); i++)
	{
		if (imv[i].status() == 0)
		{
			pointnumber++;
		}
	}

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
		cout << "not enough valid image files in the directory or no valid directory ";
		if (writetofile)
		{
			myfile1 << "not enough valid image filea in the directory or no valid directory";
		}
		return -1;
	}

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
		if (binominal(pointnumber, outliers1) > 7054320)
		{
			long k = 0;
			while (binominal(pointnumber, k) <= 7054320)
			{
				k++;
			}
			outliers1 = k;
			if (outliers1 == 0)
			{
				estimator =no_rejection;
			}
		}
	}

	if (testallestimators)
	{
		if (writetofile)
		{
			cout <<"simple linear regression" << endl;
			myfile1 << endl<< endl <<"simple linear regression" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

			cout << endl << endl<< "Siegel's median slope" << endl;
			myfile1 << endl << endl<< "Siegel's median slope" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, &myfile1, true, documentationlevel, seconds, maxseconds, maxiterations, backslash);
			if (outliers1 > 0)
			{
				cout << endl << endl << "MAD estimation" << endl;
				myfile1 << endl << endl << "MAD estimation" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_decision_in_MAD_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "S estimator" << endl;
				myfile1 << endl << endl << "S estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_decision_in_S_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "Q estimator" << endl;
				myfile1 << endl << endl << "Q estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_decision_in_Q_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "T estimator" << endl;
				myfile1 << endl << endl << "T estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_decision_in_T_ESTIMATION, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "Biweight Midvariance estimator" << endl;
				myfile1 << endl << endl << "Biweight Midvariance estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_biweight_midvariance, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "Percentage based midvariance estimator" << endl;
				myfile1 << endl << endl << "Percentage based midvariance estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_percentagebased_midvariance, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "Peirce criterion" << endl;
				myfile1 << endl << endl << "Peirce criterion" << endl;
				attempt(&imv, scale, outliers1, tolerance, use_peirce_criterion, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "removes points whose squared error is above average plus 2 standard deviations" << endl;
				myfile1 << endl << endl << "removes points whose squared error is above average plus 2 standard deviations" << endl;
				attempt(&imv, scale, outliers1, 2, tolerance_multiplies_standard_deviation_of_error, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "Grubs test at 80% significance" << endl;
				myfile1 << endl << endl << "Grubs test at 80% significance" << endl;
				attempt(&imv, scale, outliers1, 0.2, tolerance_is_significance_in_Grubbs_test, &myfile1, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);
			}
		}
		else
		{
			cout  <<"simple linear regression" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

			cout <<endl<< endl<< "Siegel's median slope" << endl;
			attempt(&imv, scale, 0, tolerance, estimator, NULL, true, documentationlevel, seconds, maxseconds, maxiterations, backslash);
			if (outliers1 > 0)
			{
				cout << endl << endl << "MAD estimation" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_decision_in_MAD_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "S estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_decision_in_S_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "Q estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_decision_in_Q_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "T estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_decision_in_T_ESTIMATION, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "Biweight Midvariance estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_biweight_midvariance, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "Percentage based midvariance estimator" << endl;
				attempt(&imv, scale, outliers1, tolerance, tolerance_is_percentagebased_midvariance, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "Peirce criterion" << endl;
				attempt(&imv, scale, outliers1, tolerance, use_peirce_criterion, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "removes points whose squared error is above average plus 2 standard deviations" << endl;
				attempt(&imv, scale, outliers1, 2, tolerance_multiplies_standard_deviation_of_error, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);

				cout << endl << endl << "Grubs test at 80% significance" << endl;
				attempt(&imv, scale, outliers1, 0.2, tolerance_is_significance_in_Grubbs_test, NULL, false, documentationlevel, seconds, maxseconds, maxiterations, backslash);
			}
		}
	}
	else
	{
		if (estimator == tolerance_is_significance_in_Grubbs_test)
		{
			if (toleranceset == false)
				tolerance = 0.9;
		}
		if (estimator == tolerance_multiplies_standard_deviation_of_error)
		{
			if (toleranceset == false)
				tolerance = 2;
		}
		if (writetofile)
		{
			attempt(&imv, scale, outliers1, tolerance, estimator, &myfile1, median, documentationlevel, seconds, maxseconds, maxiterations, backslash);
		}
		else
		{
			attempt(&imv, scale, outliers1, tolerance, estimator, NULL, median, documentationlevel, seconds, maxseconds, maxiterations, backslash);
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
	return 0;
}
