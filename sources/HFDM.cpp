// Copyright(c) < 2021 > <Benjamin Schulz.

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


// HFDM is a small program is used to test a hyperbolic fit library that can be used for autofocus purposes in telescopes.
// It needs hfd data. It can either fit historical hfd data from some previous  observations, or the user can put in
// his own data. The application also asks for several curve fitting parameters. Most is self explanatory, so not much
// documentation is needed here. For the user supplied data as well as for the historical data, the application tests 
// several curve fits with various estimators. For the historical data, more and more artificial outliers are added 
// and one can observe how the different estimators reac to the presence of an increasing number of outliers.


using namespace std;
#include <iostream>
#include <vector>
#include "focusinterpolation.h"
#include <stdlib.h>     
#include <sstream>
#include <string>
#include <chrono>



void attempt2(vector<long> xv, vector<double> yv, double scale, size_t numberofpoints, size_t outliers, double tolerance, outlier_criterion estimator, double* seconds) {

	double mainslope, mainintercept;
	vector<double> returnline_x, returnline_y, removedpoints_x, removedpoints_y;
	vector<size_t>usedindices, removedindices;


	double mainerror;
	long focpos;

	auto start = std::chrono::steady_clock::now();


	if (focusposition_Regression(xv, yv, &focpos, &mainerror, &mainslope, &mainintercept, &usedindices, &returnline_x, &returnline_y, &removedindices, &removedpoints_x, &removedpoints_y, 60, 2000000, 0, scale, false, outliers, estimator, tolerance) == false)
	{
		std::cout << "focusposition2 returned  false						" << endl;
		return;
	}
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	*seconds += elapsed_seconds.count();

	std::cout << "estimated best focusposition" << endl;
	std::cout << focpos << endl;

	std::cout << "estimated error" << endl;
	std::cout << mainerror << endl;

	std::cout << "Hyperbola parameters y=sqrt(slope (x-bestfocus)^2+intercept)" << endl;
	std::cout << "intercept   ";
	std::cout << mainintercept << endl;

	std::cout << "slope  ";
	std::cout << mainslope << endl << endl;


	std::cout << "Outliers" << endl;

	std::cout << "Motorpositions       ";
	for (size_t i = 0; i < removedindices.size(); i++)
	{
		std::cout << xv[removedindices[i]];
		std::cout << ", ";
	}
	std::cout << endl;

	std::cout << "HFD                  ";
	for (size_t i = 0; i < removedindices.size(); i++)
	{
		std::cout << yv[removedindices[i]];
		std::cout << ", ";
	}
	std::cout << endl;

	std::cout << "Used points" << endl;

	std::cout << "Motorpositions       ";
	for (size_t i = 0; i < usedindices.size(); i++)
	{
		std::cout << xv[usedindices[i]] << ",	";
	}
	std::cout << endl;

	std::cout << "HFD                  ";
	for (size_t i = 0; i < usedindices.size(); i++)
	{
		std::cout << yv[usedindices[i]] << ",	";
	}
	std::cout << endl;

}



void attempt(vector<long> xv, vector<double> yv, double scale, size_t numberofpoints, size_t outliers, double tolerance, double* seconds) {


	double mainslope1, mainintercept1, mainerror1, mainslope2, mainintercept2, mainerror2;


	long focpos1, focpos2;


	auto start = std::chrono::steady_clock::now();

	if (focusposition_Regression(xv, yv, &focpos1, &mainerror1, &mainslope1, &mainintercept1, NULL, NULL, NULL, NULL, NULL, NULL, 120, 2000, 0, scale, false, 0, no_rejection, 0) == false)
	{
		std::cout << "focusposition returned false        " << endl;
		return;
	}
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	*seconds += elapsed_seconds.count();


	std::cout << "with simple linear regression " << endl;
	std::cout << "focusposition" << endl;
	std::cout << focpos1 << endl;

	std::cout << "estimated errors" << endl;
	std::cout << "from linear regression" << endl;
	std::cout << mainerror1 << endl;
	std::cout << "Hyperbola parameters y=sqrt(slope (x-bestfocus)^2+intercept)" << endl;
	std::cout << "intercept  ";
	std::cout << mainintercept1 << endl;

	std::cout << "slope  ";
	std::cout << mainslope1 << endl << endl;

	start = std::chrono::steady_clock::now();
	if (focusposition_Regression(xv, yv, &focpos2, &mainerror2, &mainslope2, &mainintercept2, NULL, NULL, NULL, NULL, NULL, NULL, 120, 2000, 0, scale, true, 0, no_rejection, 0) == false)
	{
		std::cout << "focusposition2 returned  false						" << endl;
		return;
	}
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	*seconds += elapsed_seconds.count();


	std::cout << "with robust median regression" << endl;
	std::cout << "focusposition" << endl;
	std::cout << focpos2 << endl;

	std::cout << "estimated error" << endl;

	std::cout << mainerror2 << endl;
	std::cout << "Hyperbola parameters y=sqrt(slope (x-bestfocus)^2+intercept)" << endl;
	std::cout << "b^2*d:   ";
	std::cout << mainintercept2 << endl;

	std::cout << "b^2/a^2:  ";
	std::cout << mainslope2 << endl << endl;

	std::cout << endl<< endl << "Calculation with S estimator" << endl;
	attempt2(xv, yv, scale, numberofpoints, outliers, tolerance, Least_trimmed_squares_tolerance_is_decision_in_S_ESTIMATION, seconds);


	std::cout << endl << endl << "Calculation with Q estimator" << endl;
	attempt2(xv, yv, scale, numberofpoints, outliers, tolerance, Least_trimmed_squares_tolerance_is_decision_in_Q_ESTIMATION, seconds);

	std::cout << endl << endl << "Calculation with T estimator" << endl;
	attempt2(xv, yv, scale, numberofpoints, outliers, tolerance, Least_trimmed_squares_tolerance_is_decision_in_T_ESTIMATION, seconds);

	std::cout << endl << endl << "Calculation with MAD estimator" << endl;
	attempt2(xv, yv, scale, numberofpoints, outliers, tolerance, Least_trimmed_squares_tolerance_is_decision_in_MAD_ESTIMATION, seconds);

	std::cout << endl << endl << "Calculation with Biweight-Midvarianceestimator" << endl;
	attempt2(xv, yv, scale, numberofpoints, outliers, tolerance, Least_trimmed_squares_tolerance_is_biweight_midvariance, seconds);

	std::cout << endl << endl << "Calculation with Peirce criterion " << endl;
	attempt2(xv, yv, scale, numberofpoints, outliers, tolerance, Least_trimmed_squares_use_peirce_criterion, seconds);

	std::cout << endl << endl << "Calculation with Grubbs test and significance 90 criterion " << endl;
	attempt2(xv, yv, scale, numberofpoints, outliers, 0.1, Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test, seconds);

	std::cout << endl << endl << "Calculation with standard deviation and average, outliers with 2 sigma over average get removed" << endl;
	attempt2(xv, yv, scale, numberofpoints, outliers, 2, Least_trimmed_squares_tolerance_is_significance_in_Grubbs_test, seconds);

}



int main()
{
	string input;
	cout << "Enter 1 in order to test the hyperbolic fit with your own data" << endl<< "Enter 2 or press return if you want to test the hyperbolic fit with builtin historical data (default)  ";
	cout << endl;
	std::getline(std::cin, input);
	int decision = 2;

	if (!input.empty()) {
		istringstream stream(input);
		stream >> decision;
	}

	if ((decision!=2) && (decision!=1)) {
		cout << "bad data entered";
		return -1;
	}

	if (decision == 2) {

		double tolerance = 3;
		cout << "Please enter a value for tolerance [default = 3]: ";
		string input;
		std::getline(std::cin, input);

		if (!input.empty()) {
			istringstream stream(input);
			stream >> tolerance;
		}

		cout << endl;
		cout << "Please enter a value for the maximum number of outliers that can be found (from 0 to 3 points)[default 3]";

		size_t outliers = 3;
		getline(std::cin, input);
		if (!input.empty()) {
			istringstream stream(input);
			stream >> outliers;
		}

		if ((outliers < 0) || (outliers > 3)) {
			cout << "bad data entered";
			return -1;
		}
		cout << endl;
		cout << endl;
		cout << endl;
		cout << "Please enter a value for the scale parameter. It determines the interval size where the best focus is searched." << endl << endl
			<< "scale=1 means the best focus is searched within the measured motorpositions" << endl << endl
			<< "scale=2 doubles the size of this interval and so on. default is scale=1.5 " << endl << endl
			<< "In order to see the effects of this parameter for the historical data, scale=1.5 suffices" << endl << endl
			<< "for this data, scale should be <= 10." << endl << endl;


		double scale = 1.5;
		getline(std::cin, input);
		if (!input.empty()) {
			istringstream stream(input);
			stream >> scale;
		}

		if (scale < 1) {
			cout << "bad data entered";
			return -1;
		}
		if (scale > 10) {
			cout << "bad data entered";
			return -1;
		}

		cout << endl;
		cout << "The value you entered for tolerance is " << tolerance;

		cout << endl;
		cout << "The value you entered for minimalmodelsize is " << outliers << endl << endl;

		cout << endl;
		cout << "The value you entered for scale is " << scale << endl << endl;

		double seconds = 0;

		auto start = std::chrono::steady_clock::now();

		cout << endl << endl << endl;
		cout << "example 1 " << endl;
		std::vector <long > xv1= { 285, 270, 255, 240, 225, 210, 195 };
		std::vector<double> yv1= { 7.41, 5.65, 4.31, 3.65, 4.35, 6.00, 8.27 };

		cout << "Motorpositions: 285, 270, 255, 240, 225, 210, 195" << endl;
		cout << "HFD data: 7.41, 5.65, 4.31, 3.65, 4.35, 6.00, 8.27" << endl << endl;

		attempt(xv1, yv1, scale, 7, outliers, tolerance, &seconds);

		cout << endl << endl << endl;
		cout << "example 2 " << endl;

		std::vector <long > xv2 = { 320, 310,300, 290, 280, 270, 260, 250, 240, 230, 220 };
		std::vector<double> yv2 = { 7.55, 6.31, 5.42, 4.83, 4.53, 4.33, 5.22, 6.01, 7.09, 8.29, 9.86 };

		cout << "Motorpositions: 320, 310, 300, 290, 280, 270, 260, 250, 240, 230, 220" << endl;
		cout << "HFD data: 7.55, 6.31, 5.42, 4.83, 4.53, 4.33, 5.22, 6.01, 7.09, 8.29, 9.86" << endl << endl;
		cout << "For this large dataset, we are using the RANSAC with randomly generated combinations." << endl;
		cout << "Since the dataset is larger, we allow the algorithm here to remove 5 outliers at maximum" << endl << endl;



		attempt(xv2, yv2, scale, 11, 5, tolerance, &seconds);
		cout << endl << endl << endl;
		cout << "example 3 " << endl;

		std::vector <long >	xv3 = { 315, 300, 285, 270, 255, 240, 225 };
		std::vector<double> yv3 = { 10.91, 9.33, 7.32, 5.52, 4.09, 3.48, 4.74 };

		cout << "Motorpositions: 315, 300, 285, 270, 255, 240, 225" << endl;
		cout << "HFD data: 10.91, 9.33, 7.32, 5.52, 4.09, 3.48, 4.74" << endl << endl;

		attempt(xv3, yv3, scale, 7, outliers, tolerance, &seconds);
		cout << endl << endl << endl;
		cout << "example 4 " << endl;
		std::vector <long >	xv4= { 255, 240, 225, 210, 195, 180, 165 };
		std::vector<double> yv4 = { 4.23, 3.56, 4.27, 5.99, 8.12, 10.69, 13.62 };

		cout << "Motorpositions: 255, 240, 225, 210, 195, 180, 165 " << endl;
		cout << "HFD data: 4.23, 3.56, 4.27, 5.99, 8.12, 10.69, 13.62 " << endl << endl;


		attempt(xv4, yv4, scale, 7, outliers, tolerance, &seconds);
		cout << endl << endl << endl;

		cout << "example 5 " << endl;
		std::vector <long >	 xv5 = { 285, 270, 255, 240, 225, 210, 195 };
		std::vector<double>  yv5 = { 7.61, 5.68, 4.28, 3.34, 4.07, 5.7, 7.82 };

		cout << "Motorpositions: 285, 270, 255, 240, 225, 210, 195 " << endl;
		cout << "HFD data: 7.61, 5.68, 4.28, 3.34, 4.07, 5.7, 7.82 " << endl << endl;


		attempt(xv5, yv5, scale, 7, outliers, tolerance, &seconds);
		cout << endl << endl << endl;
		cout << "example 6 " << endl;
		std::vector <long >	xv6= { 125, 140, 155, 170, 185, 200, 215 };
		std::vector<double> yv6 = { 20.02, 15.5, 13.5, 10.51, 8.47, 6.17,  4.51 };

		cout << "Motorpositions: 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data:  20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;

		cout << "focuspoint in this example is outside of the measured data. Continuation of the measurement yielded 226. " << endl;
		cout << "Use scale = 2 to see if the algorithms can find the correct point outside of the data" << endl;
		cout << "Median Regression seems to yield the best estimate" << endl;
		attempt(xv6, yv6, scale, 7, outliers, tolerance, &seconds);

		cout << endl << endl << endl;
		cout << "example 7 outlier in  data of example 6 removed by hand." << endl;
		std::vector <long >	xv7 = { 140, 155, 170, 185, 200, 215 };
		std::vector<double> yv7 = { 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: , 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data:  15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;


		attempt(xv7, yv7, scale, 6, outliers, tolerance, &seconds);
		cout << endl << endl << endl;
		cout << "example 8, more outliers added to data of example 6 added by hand. " << endl;
		cout << "It can be removed by the RANSAC algorithmand a tolerance value = 1 if the standard deviationand average method is used." << endl;
		cout << "The robust estimators S,Q, MAD and biweight Midvariance should require no adaption of the tolerance parameter." << endl;

		std::vector <long >	xv8 = { 95, 110 ,125, 140, 155, 170, 185, 200, 215 };
		std::vector<double> yv8 = { 19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;


		attempt(xv8, yv8, scale, 9, outliers, tolerance, &seconds);
		cout << endl << endl << endl;
		cout << "example 9, even more outliers added to data of example 6 added by hand and removing up to 5 outliers " << endl;
		cout << "They can be removed by a small or negative value for tolerance if it is set in sigma units if the standard deviation and average method is used." << endl;
		cout << "The other estimators should be more robust. The robust estimators S,Q, MAD and biweight Midvariance should require no adaption of the tolerance parameter." << endl;
		std::vector <long >	xv9 = { 80,95, 110 ,125, 140, 155, 170, 185, 200, 215 };
		std::vector<double> yv9 = { 19.5,19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: 80, 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.5,19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;


		attempt(xv9, yv9, scale, 10, 5, tolerance, &seconds);
		cout << endl << endl << endl;
		cout << "example 10, even more outliers added to data of example 6 added by hand and removing up to 5 outliers." << endl;
		cout << "Then one can see at which point the estimators break down." << endl;
		std::vector <long >	xv10 = { 65, 80,95, 110 ,125, 140, 155, 170, 185, 200, 215 };
		std::vector<double> yv10 = { 19.4, 19.5,19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: 65 ,80, 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.4, 19.5,19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;

		attempt(xv10, yv10, scale, 11, 5, tolerance, &seconds);
		cout << endl<< endl<< endl;

		cout << "example 11, even more outliers added to data of example 6 added by hand and removing up to 6 outliers." << endl;
		cout << "Then one can see at which point the estimators break down" << endl;
		std::vector <long >	xv11 = { 50, 65, 80,95, 110 ,125, 140, 155, 170, 185, 200, 215 };
		std::vector<double> yv11 = { 19.8, 19.4, 19.5,19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: 65 ,80, 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.4, 19.5,19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;

		attempt(xv11, yv11, scale, 12, 6, tolerance, &seconds);

		cout << endl << endl << endl;
		cout << "example 12, even more outliers added to data of example 6 added by hand and removing up to 7 outliers." << endl;
		cout << "Then one can see at which point the robust estimators break down" << endl;
		std::vector <long >	xv12 = { 35, 50, 65, 80,95, 110 ,125, 140, 155, 170, 185, 200, 215 };
		std::vector<double> yv12 = { 19.3, 19.8, 19.4, 19.5,19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: 35, 50, 65, 65 ,80, 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.3, 19.8, 19.4, 19.4, 19.5,19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;

		attempt(xv12, yv12, scale, 13, 7, tolerance, &seconds);
		
		cout << endl;
		cout << "time used" << endl;
		cout << seconds << endl;
	}
	else
	{
	double seconds = 0;
		size_t datapointnumber=0;
		std::cout << "How many data points you want to enter ";
		std::getline(cin, input);
		stringstream(input) >> datapointnumber;

		if (datapointnumber < 4) {
			cout << "bad data entered";
			return -1;
		}

		cout << endl;
		cout << "Please enter a  value for tolerance between -4 and 4: [default = 3]: ";

		double  tolerance = 3;
		std::getline(std::cin, input);

		if (!input.empty()) {
			istringstream stream(input);
			stream >> tolerance;
		}

		if ((tolerance < -4) || (tolerance > 4)) {
			cout << "bad data entered";
			return -1;
		}
		cout << endl;
		cout << "Please enter a value for the maximum number of outlier points 0 and " << datapointnumber-4<< "  [pointnumber- 4]:  " ;

		size_t outliers = 4;
		getline(std::cin, input);
		if (!input.empty()) {
			istringstream stream(input);
			stream >> outliers;
		}

		if ((outliers<0) || (outliers> datapointnumber-4)) {
			cout << "bad data entered";
			return -1;
		}
		cout << endl;
		cout << "Please enter a value for the scale parameter. It determines the interval size where the best focus is searched." << endl
			<< "scale=1 means the best focus is searched within the measured motorpositions" << endl << endl
			<< "scale=2 doubles the size of this interval and so on. [default is scale=1]." << endl << endl;
		double scale=1;
		getline(std::cin, input);
		if (!input.empty()) {
			istringstream stream(input);
			stream >> scale;
		}

		if (scale < 1) {
			cout << "bad data entered";
			return -1;
		}

		cout << endl;
		cout << endl;
		cout << "The value you entered for the datapointnumber is " << datapointnumber;
		cout << endl;
		cout << "The value you entered for minimalmodelsize is " << outliers << endl ;
		cout << "The value you entered for scale is " << scale << endl;
		cout << "The value you entered for tolerance is " << tolerance;

		cout << endl;
	

	
		std::vector<long> xv;
		std::vector<double> yv;


		long xa = 0;
		double ya = 0;

		xv.reserve(datapointnumber);
		yv.reserve(datapointnumber);

		for (size_t i = 0; i < datapointnumber; i++) {
			std::cout << "Enter motorposition " << i + 1<<"		";
			std::getline(cin, input);
			stringstream(input) >> xa;
			xv.push_back(xa);
		}
		for (size_t i = 0; i < datapointnumber; i++) {
			std::cout << "Enter hfd data " << i + 1<<"		";
			std::getline(cin, input);
			stringstream(input) >> ya;
			yv.push_back(ya);
		}
		std::cout << "the entered motorpositions are: " << endl;

		for (size_t i = 0; i < datapointnumber; i++) {

			std::cout << xv[i] << ",  ";
		}
		std::cout << endl;
		std::cout << "the hfd-data are: " << endl;
		long minfocus = xv[0], maxfocus = xv[0];
		for (size_t i = 0; i < datapointnumber; i++)
		{
			if (xv[i] > maxfocus)
			{
				maxfocus = xv[i];
			}
			if (xv[i] < minfocus)
			{
				minfocus = xv[i];
			}
		}
	


		for (size_t i = 0; i < datapointnumber; i++) {

			std::cout << yv[i] << ",  ";
		}
		std::cout << endl;

		attempt(xv, yv,scale, datapointnumber, outliers, tolerance,&seconds);
		cout << "time used" << endl;
		cout << seconds << endl;
	}




	cout << "press key to end ";
	getline(std::cin, input);
	return 1;
	
}



	

