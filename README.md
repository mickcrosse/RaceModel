# RaceModel
RaceModel is a MATLAB package for stochastic modelling of multisensory reaction times (RTs). It is suitable for analyzing empirical data or running simulations, and can handle datasets of unequal sizes and with missing values (NaNs). The toolbox can be used to build parallel models of multisensory information processing for both OR and AND task designs (e.g., the race model; Miller, 1982), as well as bisensory and trisensory paradigms. Parallel models can be generated under the assumption that RTs on separate sensory channels are stochastically independent (independent race model), perfectly negatively dependent (Miller's bound) or perfectly positively dependent (Grice's bound), and can be tested using either the vertical or horizontal method. Separate functions compute geometric measures of multisensory gain (violation), multisensory benefit (Otto et al., 2013) and modality switch effects.

RaceModel also includes a *systems factorial technology* framework for inferring system architecture and measuring the workload capacity of a system (Townsend & Nozawa, 1995). The latter can also be assessed for OR/AND tasks and bisensory/trisensory paradigms. For statistical analyses, we recommend using multivariate permutation tests with *tmax* correction. This method provides strong control of family-wise error rate, even for small sample sizes, and is much more powerful than traditional methods (Gondan, 2010). We provide a separate MATLAB toolbox for multivariate permutation testing [here](https://github.com/mickcrosse/PERMUTOOLS "PERMUTOOLS").

### Documentation
Crosse MJ, Foxe JJ, Molholm S (2019) RaceModel: A MATLAB Package for Stochastic Modelling of Multisensory Reaction Times (In prep).

## Contents
### Redundant signals (OR) task
#### Bisensory
* `ormodel.m` - compute parallel (race) model
* `orgain.m` - compute multisensory gain (violations)
* `orbenefit.m` - compute empirical and predicted benefits
* `orcapacity.m` - compute capacity coefficient and bounds
 
#### Trisensory
* `ormodel3.m` - compute parallel (race) model
* `orgain3.m` - compute multisensory gain (violations)
* `orbenefit3.m` - compute empirical and predicted benefits
* `orcapacity3.m` - compute capacity coefficient and bounds

### Exhaustive search (AND) task
#### Bisensory
* `andmodel.m` - compute parallel (AND) model
* `andgain.m` - compute multisensory gain (violations)
* `andbenefit.m` - compute empirical and predicted benefits
* `andcapacity.m` - compute capacity coefficient and bounds

#### Trisensory
* `andmodel3.m` - compute parallel (AND) model
* `andgain3.m` - compute multisensory gain (violations)
* `andbenefit3.m` - compute empirical and predicted benefits
* `andcapacity3.m` - compute capacity coefficient and bounds

### System Architecture
* `sft.m` - systems factorial technology framework
* `biasmodel.m` - compute bias model
* `biasgain.m` - compute multisensory gain (violations)
* `biasbenefit.m` - compute empirical and predicted benefits
 
### Modality Switch Effects
* `trialhistory.m` - separate RTs based on trial history
* `switchcost.m` - compute modality switch effects

### Accuracy
* `f1score.m` - compute F1 score of a test's detection accuracy
 
### Preprocessing
* `clearnrts.m` - perform outlier correction procedures
* `rt2pdf.m` - convert RTs to a probability density function
* `rt2cdf.m` - convert RTs to a cumulative distribution function
* `rt2cfp.m` - convert RTs to a cumulative frequency polygon 
* `cfp2q.m` - convert a cumulative frequency polygon to quantiles
* `getauc.m` - compute the area under the curve
