## Test environments
* local OS X install, R 3.5.2
* ubuntu 14.04 (on travis-ci), R 3.5.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

## Other notes from the Maintainer

Dear Volunteer: thank you so much for your kind work on checking the package! Please find a few notes below. 

Update notes: 

* I consulted r-package-devel mailing list and it seems the package update rejection (dated Jun 17) is a false-positive. 

* Uopdating package as introduced major speed improvement by substituting some functionality with one from `dvmisc` (c++ based). 

Technical notes: 

* Please kindly DISABLE VIGNETTE CHECKING for my package (as it has been done upon initial submission). 

* There are NO references describing the methods in the package to be included in the Description field of the DESCRIPTION file as of now. 

