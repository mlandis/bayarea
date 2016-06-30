#bayarea
#### Bayesian inference of historical biogeography for many discrete areas

This is the software implementation of the method described in our recent Systematic Biology paper:

Landis, M. J., Matzke, N. J., Moore, B. R., & Huelsenbeck, J. P. (2013). Bayesian analysis of biogeography when the number of areas is large. Systematic Biology, 62(6), 789-804. [[link]](http://sysbio.oxfordjournals.org/content/62/6/789)

#####Quick guide
Download the most recent zip file. Unzip it. Enter the unzipped folder in the command line. Just type
```
chmod +x bayarea my_run.sh ./my_run.sh
```
to run a pre-loaded example analysis.

#####Slow guide
If the included command-line program doesn't work for your computer, you will need to compile the source code to produce an executable. In general, this is an easy process -- just follow the instructions in the manual. If you have gcc installed, type
```
g++ -O3 *.cpp -o bayarea
```
from within the source code folder. If you're stumped, feel welcome to email mlandis (at) gmail (dot) com for help.

#####Source code

BayArea is tested on Mac OS X 10.8 and compiled using gcc-4.2. C++ source code is available under Downloads (stable) and Source (unstable). An executable is provided with the source code, but you may have to compile it yourself.

#####Contact

All feedback is appreciated -- share your thoughts, questions, and comments regarding the software and model by emailing mlandis (at) gmail (dot) com.

====
##Updates

**June 29, 2016**

Migrate from GoogleCode to GitHub.

This program is no longer being actively supported. Almost all functionality has been ported to [RevBayes](http://github.com/revbayes/revbayes).

v1.0.3
* Improved MCMC starting states for better chances of the chain reaching stationarity.
* Add support for tuning proposal windows for rate and distance power parameter proposals, via `-rateProposalTuner` and `-distanceProposalTuner`.
* Add support for using standard path sampling versus auxillary variable path sampling as defined in Landis et al. (2013), toggled with the `-useAuxillarySampling` flag. The standard path sampling algorithm is guaranteed to work, which roughly follows stochastic mapping algorithm proposed by Nielsen (2002), while the auxillary path sampling algorithm lacks a mathematical proof for convergence but appears to work in practice for datasets I've examined. Credit to Jan Irvahn and Vladimir Minin for identifying this issue.


**March 13, 2014**

New tutorial hosted at http://treethinkers.org/tutorials/biogeography-with-bayarea/. Give it a try!

**September 11, 2013**

(no BayArea version update)
* Numerous updates to bayarea-fig, in particular you may now assign colors to partitions of areas.

**July 9, 2013**

v1.0.2
* File handling now works for all end line symbols.
* Initialization output block now shows the random seed when left uninitialized.
* Added setting to constrain the distance power parameter to be positive (set “-geoDistancePowerPositive=T”, default False). This is useful for global and sparse distributions with high rates of dispersal, which may send the distance power to very large negative values.
* Added setting to approximate distances computations to speed up analyses for large numbers of areas (set “-geoDistanceTruncate=T”, default False).

**June 28, 2013**

v1.0.1
* Fixed issue where analyses with small numbers of areas (e.g. N<10) would randomly freeze during analysis (thanks to Julien Vien for reporting this).
* Improved input file handling.

**May 7, 2013**

v1.0.0
* Manual now available
* Default settings set to be conservative for general analysis
* Added settings for prior and proposal densities
* Added lnL scores to area_probs.txt
* Example command-line script added
* Example vireya files added

**January 8, 2013**

v0.0.1
* Initial code clean-up for reviewers
