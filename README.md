# RaceModel
RaceModel is a MATLAB-based psychophysics toolbox for testing the race model on reaction times (RTs) generated using a multisensory detection task. Bi- or tri-sensory RT distributions can be compared to the statistical facilitation predicted by the corresponding unisensory RT distributions. The race model can be generated under the assumption that RTs on separate sensory channels are statistically independent (Raab's model), perfectly negatively dependent (Miller's bound) or perfectly positively dependent (Grice's bound). 

RaceModel can perform both vertical and horizontal tests of the race model and can handle data of different sizes and with missing values. There are also several functions for computing geometric measures of multisensory gain, empirical and predicted benefits, and modality switch effects, as well as a function for computing F1 scores of detection accuracy.
 
## Functions
### Self-terminating (OR) modelling
#### Bisensory task
* racemodel.m - compute race model
* rsegain.m - compute multisensory gain
* rsebenefit.m - compute empirical and predicted benefits
 
#### Trisensory task
* racemodel3.m - compute race model
* rsegain3.m - compute multisensory gain
* rsebenefit3.m - compute empirical and predicted benefits

### Exhaustive (AND) modelling
#### Bisensory test
* waitmodel.m - compute wait model
* eseloss.m - compute multisensory loss
* esecost.m - compute multisensory cost

#### Trisensory test
* waitmodel3.m - compute wait model
* eseloss3.m - compute multisensory loss
* esecost3.m - compute multisensory cost

### Competitive (n-1) modelling
 * competemodel.m - compute benefit predicted by competive interactions
 
### Behavioral measures
* switchcost.m - compute modality switch effects
* f1score.m - compute F1 score of detection accuracy
 
### Other
* rt2cdf.m - convert reaction times to cumulative probabilities
* rt2cfp.m - convert reaction times to a cumulative frequency polygon 
* cfp2per.m - convert a cumulative frequency polygon to percentiles
* getauc.m - compute area under the curve
