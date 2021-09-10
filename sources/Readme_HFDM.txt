// HFDM.exe is a small program can be used to determine the focus point of a telescope from
// a small series of defocused images. HFDM does, however, not have an own image analysis, in contrast to FM.exe
// HFDM.exe needs half flux diameter data from a small series of images of a star field where each of the images has to be 
// defocused to a different degree.
// HFDM.exe  can either fit historical hfd data from some previous  observations, or the user can put in
// his own hfd data. The application also asks for several curve fitting parameters. Most is self explanatory, as the application
// asks for console input and explains it. So not much documentation is needed here. 
// For the user supplied data as well as for the historical data, the application tests 
// several curve fits with various estimators. For the historical data, more and more artificial outliers are added 
// and one can observe how the different estimators react to the presence of an increasing number of artificial outliers.
