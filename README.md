# RaceModel
RaceModel is a MATLAB package for stochastic modelling of multisensory reaction times (RTs). RaceModel can model RTs for both OR and AND task designs, as well as bisensory and trisensory paradigms. Models can be generated under the assumption that RTs on separate sensory channels are statistically independent, perfectly negatively dependent or perfectly positively dependent, and can be tested using either vertical or horizontal tests. RaceModel can compute geometric measures of multisensory gain, multisensory benefit and modality switch effects, and can handle variables of unequal sizes and with missing values (entered as NaNs).

### Documentation
Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB Package for Stochastic Modelling of Multisensory Reaction Times (In prep).

## Contents
### Redundant signals (OR) task
#### Bisensory
* ormodel.m - compute OR (race) model
* orgain.m - compute multisensory gain
* orbenefit.m - compute empirical and predicted benefits
* orcapacity.m - compute capacity coefficient and bounds
 
#### Trisensory
* ormodel3.m - compute OR (race) model
* orgain3.m - compute multisensory gain
* orbenefit3.m - compute empirical and predicted benefits
* orcapacity3.m - compute capacity coefficient and bounds

### Exhaustive search (AND) task
#### Bisensory
* andmodel.m - compute AND model
* andgain.m - compute multisensory gain
* andbenefit.m - compute empirical and predicted benefits
* andcapacity.m - compute capacity coefficient and bounds

#### Trisensory
* andmodel3.m - compute AND model
* andgain3.m - compute multisensory gain
* andbenefit3.m - compute empirical and predicted benefits
* andcapacity3.m - compute capacity coefficient and bounds

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
