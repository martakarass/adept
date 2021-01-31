# adept 1.0.1

* Added a `NEWS.md` file to track changes to the package.

# adept 1.1.1

* `dvmisc`package used for sliding statistics computation. 

# adept 1.2

* Replaced `future` with `parallel` package for parallel computation. 
* Added `segmentWalking()` function. It implements algorithm to segment walking stride pattern from a raw accelerometry data time-series (x,y,z) via  Adaptive Empirical Pattern Transformation (ADEPT). Default algorithm parameters are optimized for a wrist-worn sensor and were evaluated with data collected in the free-living environment.
