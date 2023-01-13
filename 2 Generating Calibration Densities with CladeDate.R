###########################################
### Calibrating Phylogenetic Trees in R ###
###########################################

### Preliminaries ###

# Set the working directory

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/CladeDate Workshop/Code")


# Install and load the CladeDate Package

library(fGarch)

library(devtools)
install_github("evolucionario/CladeDate")
library(CladeDate)


#########################################################
### Function pdate: estimation from exact fossil ages ###
#########################################################

#Given the narrow uncertainty in the age of Limnofregata, we will ignore fossil age uncertainty for now and use the midpoint of fossil time intervals in a vector.

Fossils <- c(52.0, 48.0, 29.9, 25.7, 25.7, 25.0, 16.5)


### Uniform fossil record ###

# We can then estimate the age of Suliformes using the point estimation function pdate:

pdate(Fossils)

estimate    lower    upper 
57.07143 52.00000 70.96152

#pdate returns a point estimate and 95% confidence bounds. The option KStest = TRUE can be used to test the assumption of uniformity via a Kolmogorov-Smirnov test. 

pdate(Fossils, KStest = TRUE)

Uniform distribution not rejected (Kolmogorov-Smirnov P = 0.43)
estimate    lower    upper 
57.07143 52.00000 70.96152


### Non-uniform fossil sample ###

#The "RobsonWhitlock" method does not assume sample uniformity and uses only the two oldest fossils in the point estimator tn + (tn - tn-1) (Robson & Whitlock 1964, Solow 2003), i.e. adding the last time interval duration to the age of the oldest fossil, and the approximate confidence interval is [ tn, tn + P/(1-P) (tn - tn-1) ].
 
pdate(Fossils, method = "RobsonWhitlock")

estimate    lower    upper 
      56       52      128

#The point estimate is now a bit younger, but the uncertainty is much higher (as reflected in the upper 95% limit), as expected from an estimator that uses fewer datapoints. An alternative estimator uses only the two oldest fossils but with the additional information that each fossil must come from a different descendant lineage: the “ghost-lineage” method (Norris et al. 2015). The resultant distribution is a log-logistic distribution with shape parameter 1 and scale parameter x/2 in which x is the temporal gap between these two fossils. In the Suliformes record, the second oldest fossil is Masillastega rectirostris, a fossil skull from the Messel oil shales that show clear affinities with the Suloidea (Mayr 2002) thus representing the sister group of the Fregatidae.

pdate(Fossils, method = "NorrisGhostLin")

estimate    lower    upper 
      54       52       90

#The result of adding this information is an estimate closer to the age of the oldest fossil, and reduced uncertainty. Two additional methods are just alternative parameterizations of methods already used. The "Beta" method is based on the fact that tn/  has a Beta(n,1) distribution (Wang et al. 2009) and uses the qbeta function for estimation, producing the same results as the "StraussSadler" method.  The “penultimate gap” method of Norris et al. (2015) restricts the assumption of uniformity to the two oldest fossils that are used for estimation, does not assume knowledge of subclade affinities, and produces the same results as the "RobsonWhitlock" method.

#Finally, the “optimal linear estimation” method does not assume sample uniformity yet uses more than just the two oldest fossils: it uses a weighted sum of the k oldest fossil times (Cooke 1980, Robert & Solow 2003). Optimal weights are derived from the fact that, despite the actual distribution of the entire fossil record, the distribution of the k oldest fossil times can be modelled with a Weibull distribution (Cooke 1980, Robert & Solow 2003). This method performs well under different scenarios of non-constant fossilization or recovery potential (Wang et al. 2016) and may be the method of choice for rich fossil records (n > 20) that are not distributed uniformly. 

pdate(Fossils, method = "OLE")

estimate     lower     upper 
 63.96709  52.00000 106.80009

#In the case of the Suliformes, the number of fossils in the set is already fewer than the default for k (10) so it is not surprising that the method results in high uncertainty. Another alternative to deal with non-uniform fossil times is to evaluate whether the k oldest fossils can be assumed to be uniformly distributed (using a K-S test) and then use the "StraussSadler" method on those (Cooke 1980). 



##################################################################################
### Function clade.date: obtaining empirical distributions from time intervals ###
##################################################################################

#To demonstrate this function, we will estimate the age of the eupasseres, the main subgroup of Passeriformes including suboscine perching birds (Tyranni) and the songbirds (Passeri). The oldest fossil unequivocally representing eupasseres is Wieslochia weissi, a nearly complete skeleton found in the Rauenberg clay pits in Germany (Mayr & Manegold 2006). Wieslochia is thought to belong into the suborder Tyranni due to the presence a well-developed processus procoracoideus of the coracoid and a well-developed tuberculum ligamenti collateralis ventralis of the ulna (Mayr & Manegold 2004, 2006, Claramunt & Cracraft 2015), and was recovered as a stem Tyranni in a cladistic analysis (Ksepka et al. 2019). Calcareous nannofossils and dinoflagellate cysts indicate a fossil age in the intersection of NP23 and Subzone D14a (Maxwell et al. 2016), thus between 30.2 and 32.0 Ma (Speijer et al. 2020). The first occurrence of eupasseres in the fossil record of all other continents was taken from Claramunt & Cracraft (2015). This time, we record minimum and maximum bounds of fossil ages in a matrix in which the rows are fossils, and the columns are minimum and maximum fossil age bounds, in that order. This matrix can be built by stacking the fossil ages using:

Fossils <- rbind(c(30.2, 32.0), c(23.0, 25.0), c(16, 19), c(15.5, 16.5), c(13.6, 16.0), c(11.6, 16.0), c(5.3, 11.6))

#Then, the main function is run and the result are stored in an object of class "clade.date":

Calib <- clade.date(ages = Fossils, KStest = TRUE, n = 10000)

#The execution requires only a couple of second. A summary function prints the results in a compact way:

summary.clade.date(Calib)

	Exact one-sample Kolmogorov-Smirnov test
data:  Mages
D = 0.19383, p-value = 0.9131
alternative hypothesis: two-sided

Quantiles:
   0%   50%   95% 
30.23 33.48 42.97 

Parameters of the lognormal function:
 offset meanlog   sdlog 
30.2262  1.1792  0.8587

#Because the fossil record does not depart significantly from a uniform distribution (Kolmogorov-Smirnov test p > 0.05), the use of the "StraussSadler" method (the default) is justified, but any of the other methods implemented in pdate described above can be used with clade.date by changing the method option.

#In addition to relevant quantiles, clade.date reports the estimated parameters of a standard probability function fit to the Monte Carlo sample. clade.date  fits standard probability functions commonly used in Bayesian time-tree estimation programs. Log-normal, gamma, and exponential densities used in MrBayes (Ronquist et al. 2012) and BEAST2 (Bouckaert et al. 2019) are fit using the fitdistr() function in the MASS package (Venables & Ripley 2002), whereas skew-normal and skew-student distributions used in MCMCtree (Yang 2007) are fit with specific functions in the fGarch package (Wuertz et al. 2017). In addition to requesting specific functions (i.e.  PDFfit = "lognormal") the option "best" (the default) returns the best function among log-normal, gamma, and exponential models based on the Akaike Information Criterion. In this example, the MC distribution is better fit by a log-normal density function, which estimated parameters plus an offset (30.2 Ma, the minimum age of Wieslochia or the 0% quantile) can be used to parameterize a calibration density for eupasseres in BEAST2.

#Other options in clade.date allow for specifying the quantiles to be reported (default: p = c(0, 0.5, 0.95)), the number of pseudoreplicates ( default: n = 10000 ), the output of individual replicate values (default: repvalues = TRUE), and the output of a plot (default: plot = FALSE).

#A plotting function, plots fossil times, empirical distributions, quantiles, and probability densities (Figure 3).

plot.clade.date(Calib)

#The summary function can also covert the parameters of the probability density function to the parameterization used by MrBayes:

summary.clade.date(Calib, param="mrbayes")

Parameters of the lognormal function (MrBayes format):
 offset   mean st.dev.
 65.66   69.01   4.79

#The final example illustrates the use of the optimal linear estimation method and the fit of a skew-normal distribution that can be used in MCMCtree (Figure 4):

clade.date(ages = Fossils, method="OLE", n = 10000, PDFfit="skewStudent", plot=TRUE)


