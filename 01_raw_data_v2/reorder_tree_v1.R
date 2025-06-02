# Load the package (after installation, see above).
library(ape)
library(optimx)   # optimx seems better than R's default optim()
library(GenSA)    # GenSA seems better than optimx (but slower) on 5+ parameters, 
                  # seems to sometimes fail on simple problems (2-3 parameters)
library(rexpokit)
library(cladoRcpp)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)
library(phytools)

#######################################################
# SETUP: YOUR WORKING DIRECTORY
#######################################################
# You will need to set your working directory to match your local system

# Note these very handy functions!
# Command "setwd(x)" sets your working directory
# Command "getwd()" gets your working directory and tells you what it is.
# Command "list.files()" lists the files in your working directory
# To get help on any command, use "?".  E.g., "?list.files""#0019CD"

# Set your working directory for output files
# default here is your home directory ("~")
# Change this as you like
wd = "~/GitHub/gehyra25/01_raw_data_v2/"
setwd(wd)

# Double-check your working directory with getwd()
getwd()

#######################################################
# SETUP: Extension data directory
#######################################################
# When R packages contain extra files, they are stored in the "extdata" directory 
# inside the installed package.
#
# BioGeoBEARS contains various example files and scripts in its extdata directory.
# 
# Each computer operating system might install BioGeoBEARS in a different place, 
# depending on your OS and settings. 
# 
# However, you can find the extdata directory like this:
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

# "system.file" looks in the directory of a specified package (in this case BioGeoBEARS)
# The function "np" is just a shortcut for normalizePath(), which converts the 
# path to the format appropriate for your system (e.g., Mac/Linux use "/", but 
# Windows uses "\\", if memory serves).

# Even when using your own data files, you should KEEP these commands in your 
# script, since the plot_BioGeoBEARS_results function needs a script from the 
# extdata directory to calculate the positions of "corners" on the plot. This cannot
# be made into a straight up BioGeoBEARS function because it uses C routines 
# from the package APE which do not pass R CMD check for some reason.

#######################################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE
#######################################################
# Example files are given below. To run your own data,
# make the below lines point to your own files, e.g.
# trfn = "/mydata/frogs/frogBGB/tree.newick"
# geogfn = "/mydata/frogs/frogBGB/geog.data"

#######################################################
# Phylogeny file
# Notes: 
# 1. Must be binary/bifurcating: no polytomies
# 2. No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
# 3. Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
# 4. You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the 
#    top of the tree as fossils!  This is only a good idea if you actually do have fossils in your tree,
#    as in e.g. Wood, Matzke et al. (2013), Systematic Biology.
# 5. The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of 
#    millions of years, and the tree is 1-1000 units tall. If you have a tree with a total height of
#    e.g. 0.00001, you will need to adjust e.g. the max values of d and e, or (simpler) multiply all
#    your branchlengths to get them into reasonable units.
# 6. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#######################################################
# This is the example Newick file for Hawaiian Gehyra
# "trfn" = "tree file name"
trfn = "tree.newick"

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny (plots to a PDF, which avoids issues with multiple graphics in same window):
pdffn = "tree.pdf"
pdf(file=pdffn, width=9, height=9)

tr = read.tree(trfn)
tr
plot(tr, cex=0.5)
title("Gehyra phylogeny from Paul Oliver")
axisPhylo() # plots timescale
mtext("Millions of years ago (Ma)", side=1, line=2)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


# Edit tree to look nice
tr2 = tr
tr2 = phytools::rotateNodes(tree=tr2, nodes="all")
tr2 = read.tree(file="", text=write.tree(tr2, file=""))
plot(tr2, cex=0.5)
axisPhylo()
nodelabels(cex=0.5)
tiplabels(1:length(tr2$tip.label), cex=0.5)

tr2 = phytools::rotateNodes(tree=tr2, nodes=c(145,115))
tr2 = phytools::rotateNodes(tree=tr2, nodes=c(119))
tr2 = phytools::rotateNodes(tree=tr2, nodes=c(167))

tr2 = read.tree(file="", text=write.tree(tr2, file=""))
plot(tr2, cex=0.5)
axisPhylo()
nodelabels(cex=0.5)
tiplabels(1:length(tr2$tip.label), cex=0.5)

write.tree(tr2, file="tree2_reordered_for_viewing.newick")
write.tree(tr2, file="tree2.newick")

dev.off()


pdffn = "tree2.pdf"
pdf(file=pdffn, width=9, height=9)

tr2 = read.tree("tree2.newick")
tr2
plot(tr2, cex=0.5)
title("Gehyra phylogeny from Paul Oliver (reordered)")
axisPhylo() # plots timescale
mtext("Millions of years ago (Ma)", side=1, line=2)

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)
