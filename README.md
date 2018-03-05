# PBFE
*Predictability Bodysize Fluctuation Experiment*


## Current rationale of the Analysis (05.03.2018)

### DATA SELECTION

* Spirostomum sp. data is discarded, sporadic abundance does not allow for EDM attractor reconstruction
* also discarding incomplete dissolved Carbon and Nitrogen measurements
* summarising cryptic flagellate data (particles observed by FlowCAM together, formerly green-white and small-white particles)
* due to small time series: appending time series from different jars within the same incubator together (STITCHED.TS) (Hsieh et al. 2008), new length of the time series 198 data points (non-stitched TS referred to as INDIVIDUAL.TS)

### DATA TRANSFORMATION
*e.g CODE/DATA_PREP.R*

* based on Beninca et al. 2008, Petchey et al. 2015 and Chang et al. 2017
1. Interpolation (cubic hermite)
2. Fourth root transformation (remove sharp spikes in TS)
3. First order differencing (to obtain stationary time series, removes linear trends)
4. Rescaling (zero mean and unit variance)

(detrending not been applied because: "smoother could destroy dynamics that make signal linear" Chang 2017)

### FORECAST MODELS AND THEIR PROFICIENCY
*e.g TS_CCM, TS_PREDICT*

* use 60% of time points to train the model, predicting the rest
* use both RHO and RMSE (or nRMSE) to estimate forecast profiency
* all forecasts are based on step-wise prediction
* perform the following predictive models and compare their forecast profiencies
  * an ARIMA(0,1,0) (random walk)
  * univariate EDM 
  * multivariate EDM

(Note 1: When performing the EDM Analysis on the STITCHED.TS data the library used to train the models will be the first 40 time points of each of the composited time series, e.g if the whole time series measures (66+66+66) 198 time points, not the time points [1-160] will be taken as library but rather [1-40, 67-107, 133-173]

(Note 2: for the multivariate EDM the question emerges which variables to use when predicting a focal species. Chang et al. suggest to include variables that are directly linked to the species. Variable interactions in my foodweb are however scantly *a priori* known (mostly by vague literature) thus convergent cross-mapping was performed (CCM, Chang et al. 2017). This allows for an investigation causal relationships between two variables, which is then included to perform the multivariate EDM. E.g If variable A is influenced by variable B,D and F, as shown by CCM, these three variables will also be used in order to perform the multivariate EDM prediction for A (seperately for each replicate/incubator). 

### BODYSIZE (GROWTHRATE?) VS. FORECAST PROFIENCY

