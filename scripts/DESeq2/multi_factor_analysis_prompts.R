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

## install packages, if necessary

# import libraries


##
# Data Setup
##

# import gene count data

# check out the number of imported genes

# import grouping factors

# verify that the order of the samples in the counts and groupings files match

# convert the grouping data into factors 

# verify that the data is now of the factor type

# create DESeqDataSet list object

# inspect the list object

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
# Plotting Palettes
##

# change the graphical parameters

# view all available ghibli palettes

# close the plot and return the display to the default graphical parameters

# retrieve the vector of colors associated with PonyoMedium

# retrieve the vector of colors associated with YesterdayMedium

# retrieve the vector of colors associated with KikiMedium

##
# Data Exploration
##

# vst the data

# visualize the overall effect of the experimental conditions or any batch effects

# save the PCA

# retrieve the first two PC's % variance

# customize the PCA

# store the PCA plot

# save the PCA plot

# transpose of the transformed count matrix to get sample-to-sample distances

# convert the distances to a matrix

# update the column names

# create a list of continuous colors

# view the similarities and dissimilarities between samples

# store the clustering plot as a ggplot object

# save the plot to a png file


##
# DE Analysis Contrasts
##

# standard differential expression analysis steps are wrapped into a single function

# extract a results table with log2 fold changes, p values and adjusted p values

# check out the results

# directly specify the comparison

# check out the results

# summarize some basic tallies

# retrieve information about which variables and tests were used

# order our results table by the smallest p value

# save the ordered results to a csv file

# check the number of adjusted p-values were less than 0.5

# set the adjusted p-value cut off to 0.5

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

# show the log2 fold changes 

# save the plot to a jpg file

# shrink the log2 fold changes to remove the noise 

# inspect the shrunken log2 fold changes

# it is more useful to visualize the MA-plot for the shrunken log2 fold changes

# save the plot to a jpg file

# it can also be useful to examine the counts of reads for a single gene across the groups

# save the plot of counts for a single gene across the groups

# customize the plot of counts for a single gene

# store the plot of counts

# save the plot to a png file
