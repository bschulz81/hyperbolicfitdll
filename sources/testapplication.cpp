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


// This small program is used to test a hyperbolic fit library that can be used for autofocus purposes in telescopes

using namespace std;
#include <iostream>
#include <vector>
#include "hyperbolicfitdll.h"
#include <stdlib.h>     
#include <sstream>
#include <string>
#include <chrono>





void attempt(vector<long> xv, vector<double> yv, double scale, size_t numberofpoints, size_t outliers, double tolerance) {


	double mainslope2,mainintercept2;
	vector<double> returnline_x1, returnline_y1, removedpoints_x1, removedpoints_y1, returnline_x2, removedpoints_x2, returnline_y2, removedpoints_y2;
	vector<size_t>usedindices1, usedindices2, removedindices1, removedindices2, removedindices3, removedindices3a, removedindices4, removedindices5, removedindices6, removedindices7;


	double mainerror,mainerror2,mainerror3,mainerror4,mainerror5, mainerror5a, mainerror6,mainerror7,mainerror8,mainerror9;
	long focpos, focpos2,focpos3,focpos4,focpos5, focpos5a, focpos6,focpos7,focpos8,focpos9;


	if (focusposition_Regression(xv, yv, &focpos,&mainerror,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,120,2000,0,scale,false,0, no_rejection,0) == false)
	{
		std::cout << "focusposition returned false        " << endl;
		return;
	}
	std::cout << "estimated best focuspositions" << endl;
	std::cout << "with simple linear regression " << endl;
	std::cout << focpos << endl;


	if (focusposition_Regression(xv, yv, &focpos2,&mainerror2, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 120, 2000,0,scale,true,0, no_rejection,0) == false)
	{
	std::cout << "focusposition2 returned  false						" << endl;
	return;
	}
	std::cout << "with robust median regression"<< endl;
	std::cout << focpos2 << endl;
	
	if (focusposition_Regression(xv, yv, &focpos3, &mainerror3, &mainslope2, &mainintercept2,&usedindices1, &returnline_x1, &returnline_y1, &removedindices1, &removedpoints_x1, &removedpoints_y1,60, 2000000,0,scale,false,outliers,tolerance_is_decision_in_MAD_ESTIMATION,tolerance) == false)
	{
		std::cout << "focusposition2 returned  false						" << endl;
		return;
	}
	std::cout << "with ransac and MAD estimator"<<endl;
	std::cout << focpos3 << endl;

	if (focusposition_Regression(xv, yv, &focpos9, &mainerror9, NULL, NULL, NULL, NULL, NULL, &removedindices7, NULL, NULL, 60, 2000000, 0, scale, false, outliers, tolerance_is_biweight_midvariance, tolerance) == false)
	{
		std::cout << "focusposition2 returned  false						" << endl;
		return;
	}
	std::cout << "with biweight Midvariance" << endl;
	std::cout << focpos9 << endl;


	if (focusposition_Regression(xv, yv, &focpos4, &mainerror4, NULL, NULL, NULL, NULL, NULL, &removedindices2, NULL, NULL, 60, 2000000, 0, scale, false, outliers, tolerance_is_decision_in_Q_ESTIMATION, tolerance) == false)
	{
		std::cout << "focusposition2 returned  false						" << endl;
		return;
	}
	std::cout << "with ransac and Q estimator" << endl;
	std::cout << focpos4 << endl;

	if (focusposition_Regression(xv, yv, &focpos5, &mainerror5, NULL, NULL, NULL, NULL, NULL, &removedindices3, NULL, NULL, 60, 2000000, 0, scale, false, outliers, tolerance_is_decision_in_S_ESTIMATION, tolerance) == false)
	{
		std::cout << "focusposition2 returned  false						" << endl;
		return;
	}
	std::cout << "with ransac and S estimator" << endl;
	std::cout << focpos5 << endl;


	if (focusposition_Regression(xv, yv, &focpos5a, &mainerror5a, NULL, NULL, NULL, NULL, NULL, &removedindices3a, NULL, NULL, 60, 2000000, 0, scale, false, outliers, tolerance_is_decision_in_T_ESTIMATION, tolerance) == false)
	{
		std::cout << "focusposition2 returned  false						" << endl;
		return;
	}
	std::cout << "with ransac and T estimator" << endl;
	std::cout << focpos5a << endl;


	if (focusposition_Regression(xv, yv, &focpos6, &mainerror6, NULL,NULL, NULL, NULL, NULL, &removedindices4, NULL, NULL , 60, 2000000, 0,scale,false,outliers,tolerance_multiplies_standard_deviation_of_error,tolerance) == false)
	{
		std::cout << "focusposition2 returned  false						" << endl;
		return;
	}
	std::cout << "with ransac and average and standard deviation"<< endl;
	std::cout << focpos6 << endl;


	if (focusposition_Regression(xv, yv, &focpos7, &mainerror7, NULL, NULL, NULL, NULL, NULL, &removedindices5,NULL,NULL, 60, 2000000, 0, scale, false, outliers, use_peirce_criterion, tolerance) == false)
	{
		std::cout << "focusposition2 returned  false						" << endl;
		return;
	}
	std::cout << "with Peirce criterion" << endl;
	
	std::cout << focpos7 << endl;

	if (focusposition_Regression(xv, yv, &focpos8, &mainerror8, NULL, NULL, NULL, NULL, NULL, &removedindices6, NULL, NULL, 60, 2000000, 0, scale, false, outliers, tolerance_is_significance_in_Grubbs_test, 0.2) == false)
	{
		std::cout << "focusposition2 returned  false						" << endl;
		return;
	}
	std::cout << "with Grubb's test at 0.1 significance" << endl;
	std::cout << focpos8 << endl<<endl<<endl<<endl<<endl;
		




	std::cout << "estimated errors" << endl;

	std::cout << "from linear regression"<< endl;
	std::cout << mainerror << endl;

	std::cout << "from robust median regression"<< endl;
	std::cout << mainerror2 << endl;

	std::cout << "with ransac and MAD estimator" << endl;
	std::cout << mainerror3 <<"  Number of outliers: " << removedindices1.size()<< endl;

	std::cout << "from ransac biweight Midvariance" << endl;
	std::cout << mainerror9 << "  Number of outliers: " << removedindices7.size() << endl;

	std::cout << "with ransac and Q estimator" << endl;
	std::cout << mainerror4 << "  Number of outliers: " << removedindices2.size() << endl;

	std::cout << "with ransac and S estimator" << endl;
	std::cout << mainerror5<< "  Number of outliers: " << removedindices3.size() << endl;

	std::cout << "with ransac and T estimator" << endl;
	std::cout << mainerror5a << "  Number of outliers: " << removedindices3a.size() << endl;

	std::cout << "with ransac and average and standard deviation" << endl;
	std::cout << mainerror6 << "  Number of outliers: " << removedindices4.size() << endl;

	std::cout << "from ransac and Peirce criterion" << endl;
	std::cout << mainerror7<<"  Number of outliers: " << removedindices5.size() << endl;

	std::cout << "from ransac and Grubb's test at 0.1 significance" << endl;
	std::cout << mainerror8 << "  Number of outliers: " << removedindices6.size() << endl;






	cout << endl << endl << endl;
	cout << "example data for MAD estimator" << endl;

	std::cout << "Hyperbola parameters y=sqrt(slope (x-bestfocus)^2+intercept)" << endl;
	std::cout << "b^2*d:   ";
	std::cout << mainintercept2 << endl;

	std::cout << "b^2/a^2:  ";
	std::cout << mainslope2 << endl << endl;


	std::cout << "Outliers in original form" << endl;

	std::cout << "Motorpositions       " ;
	for (size_t i = 0; i < removedindices1.size(); i++)
	{
		std::cout << xv[removedindices1[i]];
		std::cout << ", ";
	}
	std::cout << endl;

	std::cout << "HFD                  "  ;
	for (size_t i = 0; i < removedindices1.size(); i++)
	{
		std::cout << yv[removedindices1[i]];
		std::cout << ", ";
	}
	std::cout << endl;


	std::cout << "Outliers in line coordinates"	<<endl;

	std::cout << "Motorpositions       " ;
	for (size_t i = 0; i < removedindices1.size(); i++)
	{
		std::cout << removedpoints_x1[i] << ",	";
	}
	std::cout << endl;

	std::cout << "HFD	            ";
	for (size_t i = 0; i < removedindices1.size(); i++)
	{
		std::cout << removedpoints_y1[i] << ",	";
	}
	std::cout << endl;

	std::cout << "Used points in original form" << endl;

	std::cout << "Motorpositions       ";
	for (size_t i = 0; i < usedindices1.size(); i++)
	{
		std::cout << xv[usedindices1[i]]	<< ",	";
	}
	std::cout << endl;

	std::cout << "HFD                  ";
	for (size_t i = 0; i < usedindices1.size(); i++)
	{
		std::cout << yv[usedindices1[i]] << ",	";
	}
	std::cout << endl;


	std::cout << "Used points in line coordinates" << endl;

	std::cout << "Motorpositions       ";
	for (size_t i = 0; i < usedindices1.size(); i++)
	{
		std::cout << returnline_x1[i] << ",	";
	}
	std::cout << endl;

	std::cout << "HFD	            ";
	for (size_t i = 0; i < usedindices1.size(); i++)
	{
		std::cout << returnline_y1[i] << ",	";
	}
	std::cout << endl<<endl<<endl;
	
}
inline void copytovector(long x[], double y[], size_t datapointnumnber, vector<long>* xv, vector<double>* yv) {
	(*xv).clear();
	(*yv).clear();
	(*xv).resize(datapointnumnber);
	(*yv).resize(datapointnumnber);

	for (size_t t = 0; t < datapointnumnber; t++) {
		(*xv)[t]=(x[t]);
		(*yv)[t]=(y[t]);
	}
}


int main()
{


	std::vector <long > xv;
	std::vector<double> yv;
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

		double tolerance = 2;
		cout << "Please enter a value for tolerance [default = 2]: ";
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

		if ((outliers<0) || (outliers>3)) {
			cout << "bad data entered";
			return -1;
		}
		cout << endl;
		cout << endl;
		cout << endl;
		cout << "Please enter a value for the scale parameter. It determines the interval size where the best focus is searched." << endl << endl
			<< "scale=1 means the best focus is searched within the measured motorpositions" << endl << endl
			<< "scale=2 doubles the size of this interval and so on. default is scale=1.5 "<< endl<< endl
			<< "In order to see the effects of this parameter for the historical data, scale=1.5 suffices" << endl << endl
			<< "for this data, scale should be <= 10." << endl << endl;

		
		double scale=1.5;
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
		size_t counter1 = 0;
		auto start = std::chrono::steady_clock::now();


		cout << "example 1 " << endl;
		long x[] = { 285, 270, 255, 240, 225, 210, 195 };
		double y[] = { 7.41, 5.65, 4.31, 3.65, 4.35, 6.00, 8.27 };

		cout << "Motorpositions: 285, 270, 255, 240, 225, 210, 195" << endl;
		cout << "HFD data: 7.41, 5.65, 4.31, 3.65, 4.35, 6.00, 8.27" << endl << endl;

		copytovector(x, y, 7, &xv, &yv);
		attempt(xv, yv,scale, 7, outliers, tolerance);

		cout << endl;

		cout << "example 2 " << endl;

		long x1[] = { 320, 310, 300, 290, 280, 270, 260, 250, 240, 230, 220 };
		double y1[] = { 7.55, 6.31, 5.42, 4.83, 4.53, 4.33, 5.22, 6.01, 7.09, 8.29, 9.86 };

		cout << "Motorpositions: 320, 310, 300, 290, 280, 270, 260, 250, 240, 230, 220" << endl;
		cout << "HFD data: 7.55, 6.31, 5.42, 4.83, 4.53, 4.33, 5.22, 6.01, 7.09, 8.29, 9.86" << endl << endl;
		cout << "For this large dataset, we are using the RANSAC with randomly generated combinations." << endl;
	    cout << "Since the dataset is larger, we allow the algorithm here to remove 5 outliers at maximum"<< endl << endl;

		copytovector(x1, y1, 11, &xv, &yv);
		attempt(xv, yv,scale, 11, 5, tolerance);
		cout << endl;
		
		cout << "example 3 " << endl;

		long x2[] = { 315, 300, 285, 270, 255, 240, 225 };
		double y2[] = { 10.91, 9.33, 7.32, 5.52, 4.09, 3.48, 4.74 };

		cout << "Motorpositions: 315, 300, 285, 270, 255, 240, 225" << endl;
		cout << "HFD data: 10.91, 9.33, 7.32, 5.52, 4.09, 3.48, 4.74" << endl << endl;

		copytovector(x2, y2, 7, &xv, &yv);
		attempt(xv, yv,scale, 7, outliers, tolerance);
		cout << endl;
		
		cout << "example 4 " << endl;
		long x3[] = { 255, 240, 225, 210, 195, 180, 165 };
		double y3[] = { 4.23, 3.56, 4.27, 5.99, 8.12, 10.69, 13.62 };

		cout << "Motorpositions: 255, 240, 225, 210, 195, 180, 165 " << endl;
		cout << "HFD data: 4.23, 3.56, 4.27, 5.99, 8.12, 10.69, 13.62 " << endl << endl;

		copytovector(x3, y3, 7, &xv, &yv);
		attempt(xv, yv,scale, 7, outliers, tolerance);
		cout << endl;

		
		cout << "example 5 " << endl;
		long x4[] = { 285, 270, 255, 240, 225, 210, 195 };
		double y4[] = { 7.61, 5.68, 4.28, 3.34, 4.07, 5.7, 7.82 };

		cout << "Motorpositions: 285, 270, 255, 240, 225, 210, 195 " << endl;
		cout << "HFD data: 7.61, 5.68, 4.28, 3.34, 4.07, 5.7, 7.82 " << endl << endl;


		copytovector(x4, y4, 7, &xv, &yv);
		attempt(xv, yv,scale, 7, outliers, tolerance);
		cout << endl;

		cout << "example 6 " << endl;
		long x5[] = { 125, 140, 155, 170, 185, 200, 215 };
		double y5[] = { 20.02, 15.5, 13.5, 10.51, 8.47, 6.17,  4.51 };

		cout << "Motorpositions: 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data:  20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;

		cout << "focuspoint in this example is outside of the measured data. Continuation of the measurement yielded 226. " << endl;
		cout << "Use scale = 2 to see if the algorithms can find the correct point outside of the data" << endl;
		cout << "Median Regression seems to yield the best estimate" << endl;
		copytovector(x5, y5, 7, &xv, &yv);
		attempt(xv, yv,scale, 7, outliers, tolerance);


		cout << "example 7 outlier in  data of example 6 removed by hand." << endl;
		long x5a[] = { 140, 155, 170, 185, 200, 215 };
		double y5a[] = { 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: , 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data:  15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;


		copytovector(x5a, y5a, 6, &xv, &yv);
		attempt(xv, yv,scale, 6, outliers, tolerance);
		cout << endl;

		cout << "example 8, more outliers added to data of example 6 added by hand. " << endl;
		cout << "It can be removed by the RANSAC algorithmand a tolerance value = 1 if the standard deviationand average method is used." << endl;
		cout<< "The robust estimators S,Q, MAD and biweight Midvariance should require no adaption of the tolerance parameter." << endl;

		long x5b[] = {95, 110 ,125, 140, 155, 170, 185, 200, 215 };
		double y5b[] = {19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;


		copytovector(x5b, y5b, 9, &xv, &yv);
		attempt(xv, yv, scale, 9, outliers, tolerance);

		cout << "example 9, even more outliers added to data of example 6 added by hand and removing up to 5 outliers " << endl;
		cout << "They can be removed by a small or negative value for tolerance if it is set in sigma units if the standard deviation and average method is used." << endl;
		cout<< "The other estimators should be more robust. The robust estimators S,Q, MAD and biweight Midvariance should require no adaption of the tolerance parameter." << endl;
		long x5c[] = { 80,95, 110 ,125, 140, 155, 170, 185, 200, 215 };
		double y5c[] = { 19.5,19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: 80, 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.5,19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;

		copytovector(x5c, y5c, 10, &xv, &yv);
		attempt(xv, yv, scale, 10, 5, tolerance);

		cout << "example 10, even more outliers added to data of example 6 added by hand and removing up to 5 outliers."<<endl;
		cout << "Then one can see at which point the estimators break down." << endl;
		long x5d[] = { 65, 80,95, 110 ,125, 140, 155, 170, 185, 200, 215 };
		double y5d[] = { 19.4, 19.5,19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: 65 ,80, 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.4, 19.5,19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;

		copytovector(x5d, y5d, 11, &xv, &yv);
		attempt(xv, yv, scale, 11, 5, tolerance);

		cout << "example 11, even more outliers added to data of example 6 added by hand and removing up to 6 outliers." << endl;
		cout<< "Then one can see at which point the estimators break down" << endl;
		long x5e[] = {50, 65, 80,95, 110 ,125, 140, 155, 170, 185, 200, 215 };
		double y5e[] = {19.8, 19.4, 19.5,19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: 65 ,80, 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.4, 19.5,19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;

		copytovector(x5e, y5e, 12, &xv, &yv);
		attempt(xv, yv, scale, 12, 6, tolerance);


		cout << "example 12, even more outliers added to data of example 6 added by hand and removing up to 7 outliers." << endl;
		cout <<"Then one can see at which point the robust estimators break down" << endl;
		long x5f[] = {35, 50, 65, 80,95, 110 ,125, 140, 155, 170, 185, 200, 215 };
		double y5f[] = {19.3, 19.8, 19.4, 19.5,19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51 };

		cout << "Motorpositions: 35, 50, 65, 65 ,80, 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.3, 19.8, 19.4, 19.4, 19.5,19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;

		copytovector(x5f, y5f, 13, &xv, &yv);
		attempt(xv, yv, scale, 13, 7, tolerance);

		cout << "example 13, even more outliers added to data of example 6 added by hand and removing up to 11 outliers. " << endl;
		cout <<"Then one can see at which point even the robust estimators break down " << endl;
		long x5f1[] = { 35, 50, 65, 80,95, 110 ,125, 140, 155, 170, 185, 200, 215,230,245,260,275,290,305,320,335,350,365,380 };
		double y5f1[] = { 19.3, 19.8, 19.4, 19.5,19.8,20.7 ,20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51,5.12,7.81,9.82,9.87,9.83,9.89,9.6,9.8,9.2,9.5,9.6 };

		cout << "Motorpositions: 35, 50, 65, 65 ,80, 95, 110, 125, 140, 155, 170, 185, 200, 215 " << endl;
		cout << "HFD data: 19.3, 19.8, 19.4, 19.4, 19.5,19.8,20.7, 20.02, 15.5, 13.5, 10.51, 8.47, 6.17, 4.51" << endl << endl;

		copytovector(x5f1, y5f1, 24, &xv, &yv);
		attempt(xv, yv, scale, 24, 11, tolerance);

		cout << endl;
		
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		seconds = elapsed_seconds.count();

		cout << "time used" << endl;
		cout << seconds<< endl;

	}
	else
	{

		size_t datapointnumber=0;
		std::cout << "How many data points you want to enter ";
		std::getline(cin, input);
		stringstream(input) >> datapointnumber;

		if (datapointnumber < 4) {
			cout << "bad data entered";
			return -1;
		}

		cout << endl;
		cout << "Please enter a  value for tolerance between -4 and 4: [default = 1]: ";

		double  tolerance = 1;
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
		cout << "Please enter a value for the maximum number of outlier points 0 and " << datapointnumber-4<< "  [default 4]:  " ;

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
			<< "scale=2 doubles the size of this interval and so on. [default is scale=1]." << endl << endl
			<< "If min1 is the minimum motor position and max1 is the maximum motor position." << endl
			<< "then let  middle = (min1 + max1) / 2" << endl
			<< "and min2 = middle - (middle - min1) * | scale |" << endl
			<< "and max2 = middle + (max1 - middle) * | scale | . " << endl << endl;
		cout << "currently, min1, max1 and scale should be set such that max2-min2 <= 46340. So set scale appropriately" << endl;

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
	

	

		xv.clear();

		yv.clear();

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

		attempt(xv, yv,scale, datapointnumber, outliers, tolerance);
	}

	


	cout << "press key to end ";
	getline(std::cin, input);
	return 1;
}



	

