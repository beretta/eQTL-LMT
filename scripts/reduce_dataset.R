# Package used for data sampling
library(unbalanced)
# Load the data from the csv file
args <- commandArgs(trailingOnly = TRUE)
dataset <- paste(args[1], "csv", sep = ".")
data <- read.csv(dataset, sep = ',', header = TRUE)
# Read the number of column 
n <- ncol(data)
# Remove Instances with Missing Data
data <- na.omit(data)
# Store in output the target column 
output <- data$Truth
# Store in input the independent variables' columns
input <- data[ ,-n]
# Call the OSS function. One Side Selection is an undersampling method resulting from the application of Tomek links
# followed by the application of Condensed Nearest Neighbor.
# Arguments
# X the input variables of the unbalanced dataset.
# Y the response variable of the unbalanced dataset. 
# IMPORTANT: the response variable must be a binary factor where the majority class is coded as 0 and the minority as 1.
dataReduced <- ubOSS(X = input, Y = output)
# Remove additional instances to achieve a perfectly balanced dataset
data_balanced <- ubUnder(dataReduced$X, dataReduced$Y, perc = 50, method = "percPos", w = NULL)
# Transform the previous output into a data frame
reduced <- data.frame(data_balanced$X, Class = data_balanced$Y)
# This is a "security check": print some statistics about the class variable
summary(reduced$Class)
# Export the balanced dataset into a CSV file without header
outfile <- paste(args[1], "reduced.csv", sep = "_")
write.table(reduced, outfile, row.names=FALSE, sep=",")