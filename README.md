# RaceModel
RaceModel is a MATLAB package for stochastic modelling of multisensory reaction times (RTs). RaceModel can model RTs for both OR and AND task designs, as well as bisensory and trisensory paradigms. Models can be generated under the assumption that RTs on separate sensory channels are statistically independent, perfectly negatively dependent or perfectly positively dependent, and can be tested using either vertical or horizontal tests. RaceModel can compute geometric measures of multisensory gain, multisensory benefit and modality switch effects, and can handle variables of unequal sizes and with missing values (entered as NaNs).

### Documentation
Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB Package for Stochastic Modelling of Multisensory Reaction Times (In prep).

## Contents
### Redundant signals (OR) task
#### Bisensory
* ormodel.m - compute OR (race) model
* orgain.m - compute multisensory gain (violations)
* orbenefit.m - compute empirical and predicted benefits
* orcapacity.m - compute capacity coefficient and bounds
 
#### Trisensory
* ormodel3.m - compute OR (race) model
* orgain3.m - compute multisensory gain (violations)
* orbenefit3.m - compute empirical and predicted benefits
* orcapacity3.m - compute capacity coefficient and bounds

### Exhaustive search (AND) task
#### Bisensory
* andmodel.m - compute AND model
* andgain.m - compute multisensory gain (violations)
* andbenefit.m - compute empirical and predicted benefits
* andcapacity.m - compute capacity coefficient and bounds

#### Trisensory
* andmodel3.m - compute AND model
* andgain3.m - compute multisensory gain (violations)
* andbenefit3.m - compute empirical and predicted benefits
* andcapacity3.m - compute capacity coefficient and bounds

### System Architecture
* sft.m - systems factorial technology framework
* biasmodel.m - compute bias model
* biasgain.m - compute multisensory gain (violations)
* biasbenefit.m - compute empirical and predicted benefits
 
### Modality Switch Effects
* trialhistory.m - separate RTs based on trial history
* switchcost.m - compute modality switch effects

### Accuracy
* f1score.m - compute F1 score of a test's detection accuracy
 
### Preprocessing
* clearnrts.m - perform outlier correction procedures
* rt2pdf.m - convert RTs to probability density function
* rt2cdf.m - convert RTs to cumulative distribution function
* rt2cfp.m - convert RTs to a cumulative frequency polygon 
* cfp2per.m - convert a cumulative frequency polygon to percentiles
* getauc.m - compute area under the curve
