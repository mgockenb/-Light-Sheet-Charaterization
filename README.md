# -Light-Sheet-Charaterization
Fit equations to the profile of a light sheet

The 'findWidthsUseDIRcmd.m' script will analyze a series of camera images and find the average FWHM of the sheet in each image. It will then plot the FWHM, as well as the scattered data set with standard deviation, of the light sheet at each position. 

The 'spliceLSimgs_fromBMPs.m' script will take a series of camera images and stitch them together by linearly interpolating between them. the result is an image of the light sheet along the optical axis.
