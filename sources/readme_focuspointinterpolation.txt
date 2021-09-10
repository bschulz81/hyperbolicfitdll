 This is an open source library that can help to autofocus telescopes. The library contains 4 functions and a class with 8 constructors and 5 functions.

One function (the earlier one) fits half flux diameter data to a hyperbola, from which the correct focus point can be interpolated. This function has no image analysis, which a third party program must provide on its own.

The new versions of the library contain an image class which can be constructed from a path to a fits file, a fits file structure or an array with image data. The image data is then analyzed in fourier space. A function is provided that takes a vector of image classes constructed from several images and fits it to a parabola. From this, the correct focus point can be determined.

The library is currently used in preview version of the software Astro Photography tool (APT) for the autofocus routine.

Two test applications are provided.

FM.exe expects a folder with fits files that are defocused at various levels as an argument along with some other parameters. From this the correct focus point is then interpolated with robust statistics and extensive information on the quality of the data is provided.

HFDM.exe has no image analysis. It can compute the focus point from hfd data of defocused images from starfields. If no own hfd data is provided, HFDM will analze HFD data from some historical observations. More and more artificial outliers are added to the data and the user can compare how the different statistical methods behave in the presence of more and more outliers.

The source code is provided in the folder /sources. A setup for windows is provided in the /binaries folder.

The library uses a RANSAC algorithm that selects various combinations of the input data for either a linear regression or a repeated median fit. In order to decide if data outside of a given combination should be added to a final configuration of points, various robust statistical outlier detection methods are used. If a point is not considered an outlier, it is added to the combination and a fit for the entire set is made. Then, the process starts from another selected combination.

The library makes use of an algorithm for student's distribution, which can be found at Smiley W. Cheng, James C. Fu, Statistics & Probability Letters 1 (1983), 223-227

The library also makes use of Peirce's outlier test. An algorithm for this method was developed in Gould, B. A, Astronomical Journal, vol. 4, iss. 83, p. 81 - 87 (1855).

The library also has the possibility to use MAD, S and Q estimators.

These estimators are extensively described in Peter J. Rousseeuw, Christophe Croux, Alternatives to the Median-Absolute Deviation J. of the Amer. Statistical Assoc. (Theory and Methods), 88 (1993),p. 1273,

Christophe Croux and Peter J.Rousseeuw, Time-effcient algorithms for two highly robust estimators of scale, In: Dodge Y., Whittaker J. (eds) Computational Statistics. Physica, Heidelberg, https :doi.org/10.1007/978-3-662-26811-7_58

The library has the option that one can use Siegel's repeated median from Siegel, Andrew (September 1980). "Technical Report No. 172, Series 2 By Department of Statistics Princeton University: Robust Regression Using Repeated Medians" within the RANSAC.

In practice, repeated median regression is a rather slow fitting method. If median regression is not used, the RANSAC will use a faster linear regression algorithm for the hyperbolic fit which was provided by Stephen King at https:aptforum.com/phpbb/viewtopic.php?p=25998#p25998).

The algorithm for image analysis in Fourier space was developed by C. Y. Tan. Its output can be fitted to a parabola from which the correct focus can be estimated.

The function focusposition_Regression2 has the following syntax:

extern "C" FOCUSINTERPOLATION_API bool focusposition_Regression2(vector* images, long* focpos, double* main_error, double* main_slope, double* main_intercept, double* theta, vector<size_t>* indices_of_used_points, vector<size_t>* indices_of_removedpoints, double stop_after_seconds=60, size_t stop_after_numberofiterations_without_improvement= 2000000, long backslash=0, double scale=1.0, bool use_median_regression=false, size_t maximum_number_of_outliers=3, outlier_criterion rejection_method= tolerance_is_decision_in_S_ESTIMATION, double tolerance=3);

It expects a pointer to a vector of image classes. From this, it interpolates the optimal focuser position from a fit with symmetric parabolas based on a RANSAC algorithm that utilizes a linear regression with a least square error comparison. If specified, the RANSAC can also use repeaded median regression. The algorithm for the linear regression that is used by the RANSAC was first published by Stephen King (username STEVE333), at https:aptforum.com/phpbb/viewtopic.php?p=25998#p25998

The function focusposition_Regression interpolates the optimal focuser position from supplied hfd data by a fit with symmetric hyperbolas based on a RANSAC algorithm that utilizes a linear regression with a least square error comparison. It functions essentially similar as focusposition_Regression2 , therefore, we first document focusposition_Regression in detail here.

focusposition_Regression has the following declaration:

extern "C" FOCUSINTERPOLATION_API bool focusposition_Regression(vector x, vector y, long* focpos, double* main_error = NULL, double* main_slope = NULL, double* main_intercept = NULL, vector<size_t>* indices_of_used_points = NULL, vector* usedpoints_line_x = NULL, vector* usedpoints_line_y = NULL, vector<size_t>* indices_of_removedpoints = NULL, vector* removedpoints_line_x = NULL, vector* removedpoints_line_y = NULL, double stop_after_seconds = 60, size_t stop_after_numberofiterations_without_improvement = 2000000, long backslash = 0, double scale = 1.5, bool use_median_regression = false, size_t maximum_number_of_outliers = 3, outlier_criterion rejection_method = tolerance_is_decision_in_S_ESTIMATION, double tolerance = 3);

it expects the following arguments:

a vector x with motor positions at which the hfd number was measured. a vector y. The hfd value of a star or a star field corresponding to motor position x [ i ] is assumed to be y [ i ] . both x and y should have at least 4 points.

If the function terminates successfully, the pointer focpos contains the estimated focusposition as a variable of type long that the function computes. This pointer must not be NULL if the function is to terminate successfully

main_error is a pointer to a double value. If the pointer is not NULL, the value of main_error will be the sum of the absolute values of the differences between the best fit and the measurement data.

main_slope and main_intercept are pointers to double values. If the pointers are not NULL, their values will be the slope and intercept of a line given by Y=slope*(x-focpos)^2+intercept. The square-root f=sqrt(Y) is the fitted hyperbola.

indices_of_used_points is a pointer to a vector which, if not NULL, will contain the indices i of the points in x [i] and y [i] that were used for the fit.

usedpoints_line_x and usedpoints_line_y are pointers to two vectors for the datatype double. If the pointers are not NULL, the vectors contain the points which are used in the fit in a coordinate system where the hyperbola is a line. One can plot them together with the line given by Y=slope*(x-focpos)^2+intercept

indices_of_removedpoints is a pointer to a vector which, if not NULL, will contain the indices i of the points in x [i] and y [i] that were not used for the fit

removedpoints_line_x and removedpoints_line_y are pointers to two vectors. If the pointers are not NULL, the vectors contain the points which were not used in the fit in the coordinate system where the hyperbola is a line. One can plot them together with the line given by Y=slope*(x-focpos)^2+intercept

The algorithm works by first selecting combination of points and fitting them to a hyperbola. This initial hyperbola is then corrected with contributions from other points. A point outside a combination is added if its error from the initial fit is deemed not to be an outlier based on various statistical methods. A new fitt with the added point is then made and the process is repeated with another initial combination of points.

The initial combination is selected randomly if the binomial coefficient of the number of points over the number of outliers is larger than 20 over 10. Otherwise, the combinations are searched deterministically.

stop_after_seconds is a parameter that stops the RANSAC after a given time in seconds has elapsed. stop_after_numberofiterations_without_improvement is a parameter that lets the RANSAC stop after it has iterated by stop_after_numberofiterations_without_improvement iterations without a further improvement of the error. Note that this parameter is not the iteration number, but it is the number of iterations without further improvement.

The parameters stop_after_seconds and stop_after_numberofiterations_without_improvement are only used if the binomial coefficient of the number of points over the number of outliers is larger than 20 over 10.

backslash is a parameter that can contain the focuser backslash in steps. The best focus position is corrected with respect to this backslash. If you already have taken account of the focuser backslash, for example by setting a suitable overshoor or a final_inwards_movement in APT or a different software or hardware correction of the backslash, set this parameter to 0

scale is a parameter of the type double that specifies the size of the interval of motor positions where the best focusposition is searched for. The default of scale is 1. Sometimes, the focus point may be outside of the interval of motor positions where the hfd was measured. let middle =(max + min) / 2 and max and min be the maximum and minimum motorposition where a hfd was measured. If an initial search finds the best focus point exactly at the right edge of the measurement interval of motor positions, then, if scale>1, the right side of the interval where the best focus is searched is enlarged. The new right boundary is then given by max = middle + (max - middle) * abs(scale) Similarly, if an initial search finds the best focus point exactly on the left side of the measurement interval, the search interval is enlarged, with the new left boundary given by min = (middle - (middle - min) * abs(scale).

use_median_regression is a parameter that specifies whether the RANSAC uses a simple linear regression or a median regression. Repeated median regression is slightly more stable against small outliers if one does not use the RANSAC algorithm.

maximum_number_of_outliers is a parameter that specifies how many outliers the ransac can maximally throw away.

rejection_method is a parameter that specifies the method which is used to reject outliers. Assume you have n datapoints. The algorithm works by searching through either all or (if the binomial coefficient of points over the number of outliers is larger than 20 over 10) randomly generated so - called minimal combinations of m=n - maximum_number_of_outliers points.

The algorithm searches for the best combination of points with the lowest error, based on linear regression, or repeated median regression. For each minimal combination, the points outside of this minimal set of m points are considered.

The error between the fit w of a minimal combination and a measurement at a motor position x is given by err_p=p(x)-w(x).

If

rejection_method==no_rejection, then the function uses every point for the fit.

rejection_method==tolerance_is_maximum_squared_error,

then a point p outside of the minimal combination are only added to the final fit if its squared error err_p*err_p fulfills

err_p*err_p<=abs(tolerance).

The RANSAC then computes an overall fit with this combination, and then constructs a new set based on a different minimal mode, until the best combination of points was found.

To specify the tolerable error directly is useful if the largest tolerable error is known, e.g. from a measurement of the seeing. By taking a series of images, an application may measure the random deviation of the hfd from the average that stars have. With outlier detection methods, one can remove the outliers of large amplitudes. setting tolerance_in_sigma_units=false, one can supply a maximally tolerable deviation from the average hfd directly to the library as the tolerance parameter. The tolerance value then corresponds to the absolute value of the maximally tolerable hfd deviation.

If rejection_method== tolerance_multiplies_standard_deviation_of_error,

then a point p outside of the best minimal combination is added to the final fit if error err_p fulfills:

(err_p-a)<= tolerance*s.

where a and s are the average and standard deviation of the errors err_p for each measured point with respect to the fit of the given minimal combination w that the algorithm tries.

If most of your data is good, then a is rather small and you may want to include most points. In that case, set tolerance=1, for example. This excludes only the largest errors that arose from seeing.

However, assume that for some motor positions far away from the focus point, the hfd values do not follow a hyperbola at all and therefore deviate much. Or assume that you have a large seeing error and you need to throw many points away.

In that case, your data has many points which have a very large error or difference when compared to a hyperbola. And only a subset of the data close to the focus point may have a small error. In order to exclude larger datasets which do not correspond to a hyperbola at all, you would only want to retain errors which are smaller than the average error a. In that case, set tolerance=-1 or even -2 if many of your points do have very large errors.

The average and standard deviation is not robust against outliers. If rejection_method== tolerance_is_standard_deviation_of_squared_error_from_mean, is used, the tolerance value would have to change depending on whether your data contains many outliers or not.

Therefore, the library provides more robust methods for accurate removal of outliers where such changes have to happen less frequently.

If

rejection_method==tolerance_is_decision_in_MAD_ESTIMATION, or rejection_method==tolerance_is_decision_in_S_ESTIMATION, or
rejection_method==tolerance_is_decision_in_Q_ESTIMATION, or rejection_method==tolerance_is_decision_in_T_ESTIMATION,

then, MAD, S, Q or T estimators are used.

A point is then added to a minimal combination if

abs(err_p-median(err_p))/ estimator <= tolerance

where median(err_p) is the median of the errors, and estimator is then the MAD, Q, S or T estimator. Q and S estimators are better suited for asymmetric distributions.

if MAD, Q, or S estimators are used. the tolerance parameter should be around 2...3.

if rejection_method==tolerance_is_significance_in_Grubbs_test,

then a point p is added to a minimal combination if its error err_p is not an outlier according to the Grubbs test. The tolerance parameter is here a significance value based on student's distribution. This means, in order to remove an outlier with 90% significance with Grubb's method, one should set tolerance=0.1, for 80% significance, one should select 0.2.

if rejection_method==use_peirce_criterion,

then a point is added to a minimal combination if it does fulfill Peirce's criterion for not being an outlier.

For this the rejection method according to Peirce, the tolerance parameter is not used. Note that number of outliers that Peirce's criterion finds is often strongly influenced by the parameter maximum_number_of_outliers

tolerance is a parameter that is used for the different rejection methods.

The function focusposition_Regression returns false on error and true if it is successful.

The function findbackslash finds the focuser backslash from two measurements of the best focus positions. It returns true if successful. The idea to correct for the backslash in this way was first published by Jim Hunt (username JST200) at https:aptforum.com/phpbb/viewtopic.php?p=26265#p26265

The function expects 2 data sets of motor positions x1 and x2 and hfd values y1 and y2. The function then calls focusposition_Regression and makes a curve fit for both sets. If the function is successful, the variable value backslash contains the focuser backslash given by the difference

optimal_focus_point_2 - optimal_focus_point1

the parameters main_slope1, main_intercept1, indicesofusedpoints1, used_points1_line_x, used_points1_line_y, indicesofremovedpoints1, removedpoints1_line_x, removedpoints1_line_y, and main_slope2, main_intercept2, indicesofusedpoints2, used_points2_line_x, used_points2_line_y, indicesofremovedpoints2, removedpoints2_line_x, removedpoints2_line_y have the same meaning as the corresponding return parameters of focusposition_Regression, just that they are for the two separate datasets x1,y1 and x2,y2.

The parameters stop_after_seconds, stop_after_numberofiterations_without_improvement, maximum_number_of_outliers, tolerance, scale, use_median_regression, use_median_regression have the same meaning as the corresponding parameters in focusposition_Regression and are used for the fits of both datasets.

The class image contains the code for image analysis. It has the following constructors:

image::image(string* filename, long focuser_position = LONG_MIN, double hfd = DBL_MIN, double fwhm = DBL_MIN); constructs the image class from a fits file given by filename. The focuser position of the image is either supplied as a parameter or read from the fits file. The class can also store hfd and fwhm values if supplied. These values are, however, not used by the curve fitting procedures. But they may be supplied by applications with an own image analysis.

image::image(fitsfile* fptr, long focuser_position, double hfd = DBL_MIN, double fwhm = DBL_MIN); constructs the image class from a fits file object. The focuser position of the image is either supplied as a parameter or read from the fits file. The class can also store hfd and fwhm values if supplied. These values are, however, not used by the curve fitting procedures. But they may be supplied by applications with an own image analysis.

image::image(size_t width, size_t height, long focuser_position, vector * imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN); constructs the image class from a double array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used by the curve fitting procedures. But they may be supplied by applications with an own image analysis.

image::image(size_t width, size_t height, long focuser_position, vector * imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN); constructs the image class from a float array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used by the curve fitting procedures. But they may be supplied by applications with an own image analysis.

image::image(size_t width, size_t height, long focuser_position, vector * imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN); constructs the image class from a long long array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used by the curve fitting procedures. But they may be supplied by applications with an own image analysis.

image::image(size_t width, size_t height, long focuser_position, vector * imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN); constructs the image class from a long array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used by the curve fitting procedures. But they may be supplied by applications with an own image analysis.

image::image(size_t width, size_t height, long focuser_position, vector * imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN); constructs the image class from a short array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used by the curve fitting procedures. But they may be supplied by applications with an own image analysis.

image::image(size_t width, size_t height, long focuser_position, vector <int8_t>* imagedata, double hfd = DBL_MIN, double fwhm = DBL_MIN); constructs the image class from an int8_t array. dimension2 is the height, dimension1 the width of the image. The focuser position of the image must be supplied. The class can also store hfd and fwhm values if supplied. These values are, however, not used by the curve fitting procedures. But they may be supplied by applications with an own image analysis.

double image::invpower();
returns the inverse power of an image in fourier mode. This is the data that is fitted to a parabola.

long image::focuser_position(); returns the focuser position of the image class

double image::hfd(); returns the hfd that can be supplied by the user in the constructor(note that this value is not used in the algorithms, but can be stored as an option, for applications which have their own hdf analysis algorithm.)

double image::fwhm(); returns the fwhm (note that this value is not used in the algorithms, but can be stored as an option, for applications which have their own hdf analysis algorithm.)

int image::status(); returns the pstatus variable. status=0 means the image class was successfully constructed.

extern "C" FOCUSINTERPOLATION_API bool focusposition_Regression2(vector* images, long* focpos, double* main_error, double* main_slope, double* main_intercept, double* theta, vector<size_t>* indices_of_used_points, vector<size_t>* indices_of_removedpoints, double stop_after_seconds, size_t stop_after_numberofiterations_without_improvement, long backslash, double scale, bool use_median_regression, size_t maximum_number_of_outliers, outlier_criterion rejection_method, double tolerance);

focusposition_Regression2 computes the focus point from a vector of validly constructed image classes which must all have status 0. these image classes can be constructed from a path to a fits file, a fits file structure, or an array with image data. The constructor then computes the power spectrum, with which focusposition_Regression2 can work. The parameters of focusposition_Regression2 are similar as in focusposition_Regression: double* main_slope, double* main_intercept are the slope and intercept for the functions invpower=slope(x-focus_point)^2+intercept, where invpower is the given by the invpower method if x is the focus motor position of an image class. ínvpower is the inverse of the power function from the fourier analysis. Currently, the parameter theta is not used. The parameters from focusposition_Regression where the fitted curve is returned in a coordinate system where it is represented by a line are omitted.

extern "C" bool FOCUSINTERPOLATION_API findbackslash_Regression(long* backslash, vector x1, vector y1, vector x2, vector y2, double* main_error1=NULL, double* main_slope1=NULL, double* main_intercept1=NULL, vector<size_t>*indicesofusedpoints1=NULL, vector*used_points1_line_x=NULL, vector*used_points1_line_y=NULL, vector<size_t>indicesofremovedpoints1=NULL, vectorremovedpoints1_line_x=NULL, vectorremovedpoints1_line_y=NULL, double main_error2=NULL, double main_slope2=NULL, double main_intercept2=NULL, vector<size_t>*indicesofusedpoints2 = NULL, vector*used_points2_line_x=NULL, vector*used_points2_line_y=NULL, vector<size_t>*indicesofremovedpoints2=NULL, vector*removedpoints2_line_x=NULL, vector*removedpoints2_line_y=NULL, double stop_after_seconds = 60, size_t stop_after_numberofiterations_without_improvement = 2000000, double scale = 1.5 ,bool use_median_regression = false, size_t maximum_number_of_outliers = 3, outlier_criterion rejection_method = tolerance_is_decision_in_S_ESTIMATION, double tolerance = 3); The function findbackslash_Regression finds the focuser backslash from two measurements of the best focus positions. It returns true if successful. The idea to correct for the backslash in this way was first published by Jim Hunt (username JST200) at https:aptforum.com/phpbb/viewtopic.php?p=26265#p26265

The function expects 2 data sets of motor positions x1 and x2 and hfd values y1 and y2. The function then calls focusposition_Regression and makes a curve fit for both sets. If the function is successful, the variable value backslash contains the focuser backslash given by the difference

optimal_focus_point_2 - optimal_focus_point1

The parameters main_slope1, main_intercept1, indicesofusedpoints1, used_points1_line_x, used_points1_line_y, indicesofremovedpoints1, removedpoints1_line_x, removedpoints1_line_y and main_slope2, main_intercept2, indicesofusedpoints2, used_points2_line_x, used_points2_line_y, indicesofremovedpoints2, removedpoints2_line_x, removedpoints2_line_y have the same meaning as the corresponding return parameters of focusposition_Regression, just that they are for the two separate datasets x1,y1 and x2,y2.

The parameters stop_after_seconds, stop_after_numberofiterations_without_improvement, maximum_number_of_outliers, tolerance, scale, use_median_regression, use_median_regression have the same meaning as the corresponding parameters in focusposition_Regression and are used for the fits of both datasets.

extern "C" FOCUSINTERPOLATION_API bool findbackslash_Regression2(long* backslash, vectorimages1, vectorimages2, double main_error1, double main_slope1, double* main_intercept1, vector<size_t>indicesofusedpoints1, vector<size_t>indicesofremovedpoints1 , double main_error2, double main_slope2, double* main_intercept, double* theta1, double* theta2, vector<size_t>*indicesofusedpoints2 , vector<size_t>*indicesofremovedpoints2 , double stop_after_seconds, size_t stop_after_numberofiterations_without_improvement, double scale , bool use_median_regression , size_t maximum_number_of_outliers , outlier_criterion rejection_method, double tolerance );

The function findbackslash_Regression2 finds the focuser backslash from two interpolations of the best focus positions. In contrast to findbackslash_Regression, it accepts 2 vectors of image classes from 2 consecutive measurements where the focus motor moved in different directions, and is otherwise smilar. The parameter theta is currently not used. The function accepts various arguments as pointers. These are filled upon return with data if the function is successful. If one does not need them, the pointers can be set to 0.

The first sample application, HFDM.exe is a command line tool that can both fits historical hfd data and hfd data which is put in by the user. Its usage is very simple as it expects user input and explains what it needs at every step. It is therefore not explained here.

Below is the documentation of the more recent test application FM.exe that can be used on real image data without needing another application.

documentation of the example program FM. exe: FM.exe is a test application for the curve fitting library that can determine the focus point of a telescope. The test application expects a folder with fits files at different focuser positions and then interpolates the focuser position where the telescope is at focus.

In order for this to work, the fits files need to have the focuser position recorded under either one of the following Keywords: FOCUSPOS, FOCUSERPOS, FOCUSERPOSITION, FOCUSPOSITION, FOCUSMOTORPOSITION, FOCUSMOTORPOS

For example a basic usage of this program may be FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w This would analyze the fits files in the given directory. It will make several curve fits for most of the estimators available in the library and write the results into a file. If one does just want to have a single curve fit with a certain estimator, e.g. the S estimator, then one may write FM.exe -d Path_To_my_fits_files -info 3 -e S -w If one just wants to print out the focus point for the fit with the S estimator and nothing else (which may be useful if this is fed into another rrogram, write FM.exe -d Path_To_my_fits_files -info 0 -e S -w the above examples use a default maximum number of outliers: It is the pointnumber (or the number of images) -7, this means that there are at least 7 points which are assumed not to be outliers. The maximum number of outliers can be specified with the -mo option, e.g FM.exe -d Path_To_my_fits_files -info 0 -e S -w -mo 3 would find at maximum 3 outliers. If this number could reduce the number of usable points below 7, it is set to 0. in general, FM.exe stops with an error if there are not at least 7 images provided. Currently FM.exe accepts the following parameters

-d Directoryname or --Directory Directoryname, where Directoryname is the path to a folder with fits files -w or --Write_to_file if the application should document its output in a file called Directoryname\results.txt -info number or --Infolevel number, where number ranges from 0 [default] to 3, which sets the level of information that the application prints. 0 means it only prints out the focus point. 3 means it prints the full information -m or --Median if the curve fit is done with Siegel's median fitting. The default is false because the median slope is somewhat slow -sc number or --Scale number, where scale is a floating point number >=1, that enlarges the area where the focus point is searched. The default is 1 -bs number or --Backslash number, where number is of type long and should be a previously measured focuser backslash that is subtracted from the estimated focus position. The default is 0. -sec number or --Seconds_for_ransac number where number is of type double and >0. It specifies the number of seconds after which the ransac is stopped. The default is 60 -it number or --Iterations_for_ransac number where number is of type long and >0. It specifies the number of iterations without improvement after which the ransac is stopped. The default is 2000000 -mo number or --Maximum_number_of_outliers number, where number is of type long>=0. It specifies the largest number of outliers that the algorithm can find. The default is the number n of supplied images minus 7. if n<7, then the number of outliers it can find is set to 0. If the binomial coefficient of n and n-8 is larger than 7054320, then the default number of outliers is given by the largest number k where the binomial coefficient of n,k is smaller than 7054320. -tol number or --Tolerance number, where number is a value of type float. It sets a cut off value for the estimators. For all estimators the default is 3.0, except for the estumators Significance_in_Grubbs_test, where it is 0.2,Standard_deviation, where it is 2.0 -e ESTIMATOR or --Estimator ESTIMATOR defines the estimator for the outlier removal procedure. ESTIMATOR can be one of the following variables: S or 5 This mmeans the S estimator is used, which is the default. The tolerance is the cut-off in |squared_error-average_squared_error|/S(squared_errors)<=tolerance ALL This mmeans that, provided the maximum number of outliers is >0) several interpolations are made where most of the estimators are tried and the different results noted. For the standard deviation, a tolderance of 2 is used, and for the Grubbs test, a tolerance value of 0.2 is used. For the other estimators, the user supplied value or the default 3 is used. If the maximum number of outliers is 0, then only simple linear regression and the median slope are tested. NO_REJECTION or 0 This means that no outliers are removed and the tolerance value is ignored. The same behavior can be set by using -mo 0 MAXIMUM_SQUARED_ERROR or 1 This means the tolerance is the maximum squared error between the fit and the measured points STANDARD_DEVIATION or 2 this means that outliers are rejected if their squared_error with respect to the fit is larger than average_squared_error+tolerance*standard_deviation_of_squared_errors SIGNIFICANCE_IN_GRUBBS_TEST or 3 This means that the tolerance value defines the significance level in the Grubbs test. MAD or 4. This means that the MAD estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/MAD(squared_errors)<=tolerance Q or 6 This means that the Q estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/Q(squared_errors)<=tolerance T or 7 This means that the T estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/T(squared_errors)<=tolerance PEIRCE_CRITERION or 8 This means that the Peirce criterion is used BIWEIGHT_MIDVARIANCE or 9 This means that the biweight midvariance estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/biweight-midvariance(squared_errors)<=tolerance PERCENTAGE_BASED_MIDVARIANCE or 10 This means that the percentage-based midvariance estimator is used and the tolerance is the cut-off in |squared_error-average_squared_error|/percentage_based_midvariance(squared_errors)<=tolerance For example a basic usage of this program may be FM.exe -d Path_To_my_fits_files -info 3 -e ALL -w This would analyze the fits files in the given directory. It will make several curve fits for most of the estimators available in the library and write the results into a file.
