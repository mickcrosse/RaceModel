# RaceModel
RaceModel is a MATLAB-based psychophysics toolbox for testing the race model on reaction times (RTs) generated using a multisensory detection task. Bi- or tri-sensory RT distributions can be compared to the statistical facilitation predicted by the corresponding unisensory RT distributions. The race model can be generated under the assumption that RTs on separate sensory channels are statistically independent (Raab's model), perfectly negatively dependent (Miller's bound) or perfectly positively dependent (Grice's bound). 

RaceModel can perform both vertical and horizontal tests of the race model and can handle data of different sizes and with missing values. There are also several functions for computing geometric measures of multisensory gain, empirical and predicted benefits, and modality switch effects, as well as a function for computing F1 scores of detection accuracy.
 
## Index
### Redundant signals (OR) task
#### Bisensory
* ormodel.m - compute OR (race) model
* orgain.m - compute multisensory gain
* orbenefit.m - compute empirical and predicted benefits
 
#### Trisensory
* ormodel3.m - compute OR (race) model
* orgain3.m - compute multisensory gain
* orbenefit3.m - compute empirical and predicted benefits

### Exhaustive search (AND) task
#### Bisensory
* andmodel.m - compute AND model
* andgain.m - compute multisensory gain
* andbenefit.m - compute empirical and predicted benefits

#### Trisensory
* andmodel3.m - compute AND model
* andgain3.m - compute multisensory gain
* andbenefit3.m - compute empirical and predicted benefits

### Bias (*n*−1) model
* biasmodel.m - compute bias model
* biasgain.m - compute multisensory gain
* biasbenefit.m - compute empirical and predicted benefits
 
### Other measures
* switchcost.m - compute modality switch effects
* f1score.m - compute F1 score (detection accuracy)
 
### Preprocessing
* rt2cdf.m - convert reaction times to cumulative probabilities
* rt2cfp.m - convert reaction times to a cumulative frequency polygon 
* cfp2per.m - convert a cumulative frequency polygon to percentiles
* getauc.m - compute area under the curve
