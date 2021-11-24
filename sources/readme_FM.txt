 Copyright(c) < 2021 > 
<Benjamin Schulz> 


Permission is hereby granted, free of charge, to any person obtaining a copy
of this softwareand associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and /or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

FM.exe is a test application for the  curve fitting library that can determine the focus point
of a telescope. The test application expects a folder with fits files at different focuser positions and then interpolates
the focuser position where the telescope is at focus. 

 A CMakeLists.txt file is provided in the sources folder. In order to compile FM and the library, 
one has to change the paths for Open-CV and libfitsio headers and libraries at the places
where it is written in the comments of the CMakeLists.txt.
Then one has to go into a console window. Change the currend directlory to the directory where the sources are found.
Finally type "cmake ." followed by "make" which should compile the application.

In order for this application to work, the fits files need to have the focuser position recorded under either one of the following 
Keywords: FOCUSPOS,FOCUSERPOS,FOCUSERPOSITION,FOCUSPOSITION,FOCUSMOTORPOSITION,FOCUSMOTORPOS

For example a basic usage of this program may be 
FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w 

This would analyze the fits files in the given directory and yield a detailed output. It will make several curve fits for most of the estimators available in the library and write the results into a file.
If one does just want to have a single curve fit with a certain estimator, e.g. the S estimator solving the least trimmed squares problem and a default tolerance, then one may write
FM.exe -d Path_To_my_fits_files -info 3 -e lts_S -w


 If one just wants to print out the focus point for the fit with the S estimator and nothing else (which may be useful if this is fed into another program, write -info 0


FM.exe -d Path_To_my_fits_files -info 0 -e lts_S 
If one wants to have more precision of the curve fit, one can specify the options -lb and -lt. This would start with an initial fit of the function

 Power=1/(alpha*(x-x0)^2+gamma)
and then use these results as initial values in a nonlinear Levenberg-Marquardt algorithm.

If both -lb and -lt are specified, the function
Power=1/(alpha*(x-x0)^2+gamma)+exp(-theta/(alpha*(x-x0)^2+gamma))/(alpha*(x-x0)^2+gamma)+beta

will be fitted. if just -lb is given, then, theta is left to zero. if just -lt is specified, then beta is left to zero.

FM.exe -d Path_To_my_fits_files -info 0 -e lts_S -lt -lb

will make an initial linear trimmed least squares fit. Then, the non-linear fit is started with the inliers from the initial fit.
The function then looks if points that were outliers in the initial fit have to be included when the non-linear fit was made. Then a final
non-linear fit is made with the resulting point set.


The above examples use a default maximum number of outliers in the least trimmed squares problem: 
It is the pointnumber (or the number of images) -7, this means that there are at least
7 points which are assumed not to be outliers. The maximum number of outliers can be specified with the -mo option, e.g

FM.exe -d Path_To_my_fits_files -info 0 -e lts_S -w -mo 3

would find at maximum 3 outliers. If this number could reduce the number of usable points below 7, it is set to 0. 
in general, FM.exe stops with an error if there are not at least 7 images provided.

the trimmed least squares problem is computationally intensive. Its complexity depends on the binomial
coefficient of the pointnumber over the maximum number of outliers.

if one has many points, the program also has the option to make a linear and non-linear regression with all points, but to throw out or linearize errors
whose values are deemed to be outliers. This is different from the trimmed least squares approach, where the entire fit is just made with subsets
of points. Here, the fit is made with all points but the error function, which is minimized, is just computed with a subset of points.

This approach can be selected with err_vanishing and err_linear estimators.

For example  
FM.exe -d Path_To_my_fits_files -info 3 -e err_linear_S  -w

would make a single linear fit with all points, but the error function it tries to minimize would consider the errors of points which
are inliers according to the S estimator as linear contributions, whereas the other errors are squared.

The command

FM.exe -d Path_To_my_fits_files -info 3 -e err_linear_S  -w -lt -lb

would increase the precision with an additional run of a nonlinear fitting algorithm for the function

Power=1/(alpha*(x-x0)^2+gamma)+exp(-theta/(alpha*(x-x0)^2+gamma))/(alpha*(x-x0)^2+gamma)+beta


Currently FM.exe accepts the following parameters

-d Directoryname  or --Directory Directoryname, where Directoryname is the path to a folder with fits files" << endl;
-w  or --Write_to_file if the application should document its output in a file called Directoryname\\results.txt" << endl;
-info number or --Infolevel number, where number ranges from 0 [default] to 3, which sets the level of information that the application prints. 0 means it only prints out the focus point. 3 means it prints the full information" << endl;
-m or --Median if the curve fit is done with Siegel's median fitting. The default is false because the median slope is somewhat slow" << endl;
-sc number or --Scale number, where scale is a floating point number >=1, that enlarges the area where the focus point is searched. The default is 1" << endl;
-bs number or --Backslash number, where number is of type long and should be a previously measured focuser backslash that is subtracted from the estimated focus position. The default is 0. " << endl;
-sec number or --Seconds_for_ransac number  where number is of type double and >0. It specifies the number of seconds after which the ransac is stopped. The default is 60" << endl;
-it number or --Iterations_for_ransac number where number is of type long and >0. It specifies the number of iterations without improvement after which the ransac is stopped. The default is 2000000" << endl;
-mo number or --Maximum_number_of_outliers number, where number is of type long>=0. It specifies the largest number of outliers that the algorithm can find if it should solve the trimmed least squares problem (the estimators with the lts_ prefix ). The default is the number n of supplied images minus 7. if n<7, then the number of outliers it can find is set to 0. If the binomial coefficient of n and n-8 is larger than 7054320, then the default number of outliers is given by the largest number k where the binomial coefficient of n,k is smaller than 7054320." << endl;
-tol number or --Tolerance number, where number is a value of type float. It sets a cut off value for the estimators. For all estimators the default is 3.0, except for the estumators Significance_in_Grubbs_test, where it is 0.2,Standard_deviation, where it is 2.0" << endl;
-lb indicates whether a Levenberg-Marquardt algorithm should be used to fit the beta parameter for the function 1/(alpha*(x-x0)^2+gamma)+beta" << endl;
-lt indicates whether a Levenberg-Marquardt algorithm should be used to fit a second order term i.e the function 1/(alpha*(x-x0)^2+gamma)+exp(-theta/(alpha*(x-x0)^2+gamma))/(alpha*(x-x0)^2+gamma) or 1 / (alpha * (x - x0) ^ 2 + gamma) + exp(-theta / (alpha * (x - x0) ^ 2 + gamma)) / (alpha * (x - x0) ^ 2 + gamma)+beta if -lb is additionally specified" << endl;

-e ESTIMATOR or --Estimator ESTIMATOR defines the estimator for the outlier removal procedure. ESTIMATOR can be one of the following variables:" << endl;
NO_REJECTION or 0 This means that no outliers are removed and the tolerance value is ignored. The same behavior can be set by using -mo 0" << endl;

ALL_lts This mmeans that the program solves the least trimmed squares problem, where the most optimal point set is searched for. Provided that the maximum number of outliers is >0) several interpolations are made where most of the estimators are triedand the different results noted.For the standard deviation,
 tolderance of 2 is used, and for the Grubbs test, a tolerance value of 0.2 is used.For the other estimators, the user supplied value or the default 3 is used.If the maximum number of outliers is 0, then only simple linear regressionand the median slope are tested.
lts_S or 5 This mmeans the least trimmed squares problem is solved where the selection criteria is the S estimator. This is  the default. The tolerance is the cut-off  in |squared_error-median_squared_error|/S(squared_errors)<=tolerance
lts_MAXIMUM_SQUARED_ERROR or 1 This means the tolerance is the maximum squared error between the fit and the measured points
lts_STANDARD_DEVIATION or 2 this means that outliers are rejected if their squared_error with respect to the fit is larger than
average_squared_error+tolerance*standard_deviation_of_squared_errors
lts_SIGNIFICANCE_IN_GRUBBS_TEST or 3 This means that the tolerance value defines the significance level in the Grubbs test. 
lts_MAD or 4. This means that the MAD estimator is used and the tolerance is the cut-off  in |squared_error-average_squared_error|/MAD(squared_errors)<=tolerance
lts_Q or 6 This means that the Q estimator is used and the tolerance is the cut-off in |squared_error-median_squared_error|/Q(squared_errors)<=tolerance
lts_T or 7 This means that the T estimator is used and the tolerance is the cut-off in |squared_error-median_squared_error|/T(squared_errors)<=tolerance
lts_PEIRCE_CRITERION or 8 This means that the Peirce criterion is used
lts_BIWEIGHT_MIDVARIANCE or 9 This means that the biweight midvariance estimator is used and the tolerance is the cut-off
in |squared_error-average_squared_error|/biweight-midvariance(squared_errors)<=tolerance;

ALL_err_vanishing This means several estimators are tested where errors (not points) are removed from the linear regression and Levenberg-Marquardt if they are deemed outliers.
err_vanishing_MAXIMUM_SQUARED_ERROR or 10 This means errors beyond a squared error given by tolerance are made to vanish
err_vanishing_STANDARD_DEVIATION or 11 this means errors are set to zero if their squared_error with respect to the fit is larger than average_abs(error)+tolerance*standard_deviation_of_abs(errors)
err_vanishing_SIGNIFICANCE_IN_GRUBBS_TEST or 12 This means that errors are set to vanish if they fail a Grubs test with a significance level given by tolerance
err_vanishing_MAD or 13. This means that errors are set to vanish if they are deemed outliers by the MAD estimator, i.e. if they fail to pass |error-median_errors|/MAD(errors)<=tolerance
err_vanishing_S or 14  This means that errors are set to vanish if they are deemed outliers by the S estimator, i.e. if they fail to pass |error - median_errors|/S(errors) <= tolerance
err_vanishing_Q or 15 This means that errors are set to vanish if they are deemed outliers by the Q estimator, i.e. if they fail to pass |error-median_errors|/Q(errors)<=tolerance
err_vanishing_T or 16 This means that errors are set to vanish if they are deemed outliers by the S estimator, i.e. if they fail to pass |error-median_errors|/T(errors)<=tolerance
err_vanishing_PEIRCE_CRITERION or 17 This means errors vanish if they do not fulfil the Peirce criterion 
err_vanishing_BIWEIGHT_MIDVARIANCE or 18 This means that errors are set to vanish if they are outliers according to the biweight-midvariance estimator is used and the tolerance is the cut-off
in |error-median_errors|/biweight-midvariance(errors)<=tolerance
ALL_err_linear This means several estimators are tested where errors (not points) are only linearly considered as abs(err) in the least squares linear regression
and Levenberg-Marquardt if they are deemed outliers.
err_linear_MAXIMUM_SQUARED_ERROR or 19 This means errors beyond a squared error given by tolerance are linearized
err_linear_STANDARD_DEVIATION or 20 this means errors are are linearized if their squared_error with respect to the fit is larger than average_abs(error)+tolerance*standard_deviation_of_abs(errors)
err_linear_SIGNIFICANCE_IN_GRUBBS_TEST or 21 This means that errors are linearized if they fail a Grubs test with a significance level given by tolerance
err_linear_MAD or 22. This means that errors are set to vanish if they are deemed outliers by the MAD estimator, i.e. if they fail to pass |error-median_errors|/MAD(errors)<=tolerance
err_linear_S or 23  This means that errors are linearized if they are deemed outliers by the S estimator, i.e. if they fail to pass |error - median_errors|/S(errors) <= tolerance
err_linear_Q or 24 This means that errors are linearized if they are deemed outliers by the Q estimator, i.e. if they fail to pass |error-median_errors|/Q(errors)<=tolerance
err_linear_T or 25 This means that errors are linearized if they are deemed outliers by the S estimator, i.e. if they fail to pass |error-median_errors|/T(errors)<=tolerance
err_linear_PEIRCE_CRITERION or 26 This means errors are linearized if they do not fulfil the Peirce criterion " 
err_linear_BIWEIGHT_MIDVARIANCE or 27 This means that errors linearized set to vanish if they are outliers according to the biweight-midvariance estimator is used and the tolerance is the cut-off in |error-median_errors|/biweight-midvariance(errors)<=tolerance
ALL means that all estimators for the ALL_lts,ALL_err_vanishing,ALL_err_linear and NO_REJECTION are tested. The curve fit for NO_REJECTION is done twice, with and without median regression for comparison.

For example a basic usage of this program may be FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w This would analyze the fits files in the given directory. It will make several curve fits for most of the estimators available in the library and write the results into a file.
