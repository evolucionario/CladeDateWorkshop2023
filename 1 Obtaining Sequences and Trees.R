###########################################
### Calibrating Phylogenetic Trees in R ###
###########################################

### Preliminaries ###

# Install Packages if Needed

#ape
#aphid

# Set the working directory

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/CladeDate Workshop/Code")



# load the ape library

library(ape)
library(phangorn)


### Obtain Sequences from Genbank ###

# Create a character vector with the Genbank Accession numbers 
AccNum <- c("AF143727","AY140770","AF143730", "AY140765","DQ137227","KT954397", "KT954398", "AF143738")

# Download sequences from GenBank
secuencias <- read.GenBank(AccNum)

# Display sequence names
names(secuencias)

# Inspect species names
attr(secuencias, "species")


# Create new names with the union of the name of the species and the Accession number and use the nnew names for the name of each sequence

names(secuencias) <- paste(attr(secuencias, "species"), names(secuencias), sep="_")

names(secuencias)



### Align Sequences ###

# load package aphid

library(aphid)

alineamiento <- align(secuencias)

alineamiento

# Save aligment in Phylip format

write.dna(alineamiento, file="Alineamiento.phy", format="sequential", nbcol=-1, indent=0, colsep="")


### Phylogenetic Tree Reconstruction using Phangorn ###

# Create an decent initial tree using the Neighbour Joining

tree <- nj(dist.dna(alineamiento))

# Create the ML model

MLfit <- pml(tree=tree, data=as.phyDat(alineamiento))

# Optimize model

MLfit <- optim.pml(MLfit, optNni=TRUE)

# Extract tree from results:

MLtree <- MLfit$tree

# Visualize tree

plot(MLtree)

# the resultant tree is unrooted so root the tree

MLtree <- root(MLtree, outgroup="Struthio_camelus_AF143727")

plot(MLtree)

# then delete the outgroup
MLtree2 <- drop.tip(MLtree, tip="Struthio_camelus_AF143727")

plot(MLtree2)


write.tree(MLtree2, file="MLtree2.tre")
