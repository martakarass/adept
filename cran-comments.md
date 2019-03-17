## Test environments
* local OS X install, R 3.5.2
* ubuntu 14.04 (on travis-ci), R 3.5.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Other notes from the Maintainer

Dear Volunteer: thank you so much for your kind work on checking the packages! Please find a few notes below. 

*  Data directory takes 4.2Mb. `acc_running` data set is mostly responsible for that size. This data contains raw subsecond-level (100 obs. per second) data collected with wearable accelerometer and as such, is, I suppose, a very rare (if not the only one) to arrive at CRAN. Please note I assume I WILL NOT BE UPDATING THIS DATA SET in future. 
    
* Doc directory 4.8Mb. The vignettes generate ggplot2 visualizations (including visualizations of high-resultion raw accelerometry data) which make the Doc directory be heavy. 

* There are NO references describing the methods in the package to be included in the Description field of the DESCRIPTION file as of now. 
