####################################################################
### Simulations Evaluating the CladeDate algorithm using Chronos ###
####################################################################


library(ape)
library(phangorn)
library(CladeDate)


# Read ML tree generated with file 2
MLtree2 <- read.tree(file="MLtree2.tre")

### Obtain a point estimate of the age of Galloanseres using CladeDate ###

## Clade A ##

Gages <- cbind(
c(65.65, 58.7, 56.4, 55.8, 51.66, 48.0, 47, 16, 0.0222),
c(66.04, 61.7, 66.04, 58.7, 51.66, 55.8, 49, 19, 0.0229))

# first estiamte using the Strauss-Sadler method

date.cladeA <- clade.date(ages= Gages, method="StraussSadler", PDFfitting=NULL, KStest=TRUE, plot=TRUE)


## Clade Suloidea ##


Sages <- matrix(c(
51.97,	51.97,
48,		48,
30.2,	29.5,
28.4,	23.03,
26,		24,
18,		15,
0.13,	0.11,
0.019,	0.001),
ncol=2, byrow=TRUE)

Sages <- Sages[,2:1]

date.cladeB <- clade.date(ages= Sages, method="StraussSadler", PDFfitting=NULL, KStest=TRUE, plot=TRUE)



### Estimate with Chronos ###

# Point calibration using the medians

# Identify the calibration node number
plot(MLtree2, cex=2); nodelabels()

# Point calibration using median values

Calib <- makeChronosCalib(phy = MLtree2,
node =c(9, 13),
age.min = c(date.cladeA$Quantiles[2], date.cladeB$Quantiles[2]) )

# Alternative using tip names to identify nodel number

Calib <- makeChronosCalib(phy = MLtree2,
node = c(getMRCA(MLtree2, tip=c("Gallus_gallus_AF143730", "Anser_albifrons_DQ137227")),getMRCA(MLtree2, tip=c("Fregata_minor_KT954397", "Sula_sula_KT954398" ))),
age.min = c(date.cladeA$Quantiles[2], date.cladeB$Quantiles[2]))


# Execute Chronos
Chrono <- chronos(phy = MLtree2, model="discrete", calibration=Calib, control=chronos.control(nb.rate.cat=5))



# To do:

# Use min-max bound instead

# Plot of chronogram and calibration densities.

# Compare point estiamtes and min-max bounds



