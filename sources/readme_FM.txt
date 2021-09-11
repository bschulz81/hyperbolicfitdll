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

// FM.exe is a test application for the  curve fitting library that can determine the focus point
// of a telescope. The test application expects a folder with fits files at different focuser positions and then interpolates
// the focuser position where the telescope is at focus. 


// A CMakeLists.txt file is provided in the sources folder. In order to compile FM and the library, 
// one has to change the paths for Open-CV and libfitsio headers and libraries at the places
// where it is written in the comments of the CMakeLists.txt.

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
