#!/usr/bin/env Rscript

##
# DE Analysis with DESeq2
##
# DESeq2 documentation
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


##
# Working Directory
##

# set the working directory


##
# Packages
##

# install packages, if necessary

# import libraries


##
# Data Setup
##
# the values in the input matrix should be un-normalized counts or  
# estimated counts of sequencing reads (SE) or fragments (PE)

# import gene count data

# check out the number of imported genes

# import grouping factors

# verify that the order of the samples in the counts and groupings files match

# convert the grouping data into factors 

# verify that the data is now of the factor type

# create DESeqDataSet list object
# the design formula expresses the variables which will be used in modeling

# inspect the list object

# by default, the reference level for factors is based on alphabetical order
# specify the reference level

# verify the re-leveling


##
# Pre-Filtering
##
# pre-filter low count genes to improve performance and improve visualizations

# here there are 3 treated samples

# keep only rows that have a count of at least 10 for a minimal number of samples

# filter the list object

# check out the number of kept genes


##
# Normalization
##
# the DESeq2 model internally corrects for library size, so transformed or 
# normalized values such as counts scaled by library size should not be used


##
# Plotting Palettes
##

# change the graphical parameters

# view all available ghibli palettes

# close the plot and return the display to the default graphical parameters

# retrieve the vector of colors associated with PonyoMedium


##
# Data Exploration
##
# for example, variance stabilizing transformations (vst) produces
# transformed data on the log2 scale which has been normalized with respect 
# to library size or other normalization factors

# vst the data

# visualize the overall effect of the experimental conditions or any batch effects

# save the PCA

# customize the PCA

# store the PCA plot

# save the PCA plot

# transpose of the transformed count matrix to get sample-to-sample distances

# convert the distances to a matrix

# update the column names

# create a list of continuous colors

# view the similarities and dissimilarities between samples
# note that we have to provide a hierarchical clustering to the heatmap  
# function based on the sample distances, or else it would be based on the 
# distances between the rows/columns of the distance matrix

# store the clustering plot as a ggplot object

# save the plot to a png file


##
# DE Analysis Contrasts
##
# with no additional arguments to results, the log2 fold change and Wald test
# p value will be for the last last level of the factored variable over the reference level

# standard differential expression analysis steps are wrapped into a single function

# extract a results table with log2 fold changes, p values and adjusted p values

# check out the results
# condition treat vs cntrl indicates that the estimates are of the 
# logarithmic fold change log2(treated/untreated)

# note that the order of the variables of the design do not matter 
# so long as the we directly specify the comparison

# check out the results

# summarize some basic tallies
# notice that the default adjusted p-value cut off is 0.1

# retrieve information about which variables and tests were used

# order our results table by the smallest p value

# save the ordered results to a csv file

# check the number of adjusted p-values were less than 0.05

# set the adjusted p-value cut off to 0.05

# summarize the results

# save the filtered results to a csv file


##
# Results Exploration
##

# retrieve the ordered row means for the top 20 most abundant genes

# setup list of colors associated with conditions

# explore the count matrix using a heatmap of the vst data

# store the clustering plot as a ggplot object

# save the plot to a png file

# show the log2 fold changes attributable to a given variable over the mean 
# of normalized counts for all the samples

# save the plot to a jpg file

# shrink the log2 fold changes to remove the noise associated with log2 fold 
# changes from low count genes without requiring arbitrary filtering thresholds

# inspect the shrunken log2 fold changes

# it is more useful to visualize the MA-plot for the shrunken log2 fold changes

# save the plot to a jpg file

# it can also be useful to examine the counts of reads for a single gene across the groups

# save the plot of counts for a single gene across the groups

# customize the plot of counts for a single gene

# store the plot of counts

# save the plot to a png file
