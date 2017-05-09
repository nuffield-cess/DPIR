

# -------------------------------------------
# -- I. Getting Started in R --
# -------------------------------------------  

# A. Basics

# Remove every objects in your working environment 
rm(list = ls())
# remove.packages()

#  Vectors 
vec1 <- rnorm(100, 75, 15)
vec2 <- as.factor(rep(c("A", "B", "C", "D"), times = 25))
vec3 <- rnorm(100, 0.5, 1)

# Vector Indices
v1[3]
v2[1:10]
v2[c(4, 8, 12)]

# 2. Matrix
m1 <- cbind(vec1, vec3)  # try cbind(vec1, vec3) 
m1  # How is this different than vec3?
class(m1)  # All columns in a matrix must be the same mode(numeric/character...), 
dim(m1)

# Referencing via Brackets : <matrix>[<row indices>,<column indices>]
m1[, 1]
m1[2, ]

# Data Frames 
mydata <- data.frame(vec1, vec2, vec3)
mydata$v1
mydata[1:10, c("v1", "v2")]

# Data Type Conversion
# as.numeric() 
# as.character()
# as.vector() 
# as.matrix() 
# as.data.frame()
# as.factor()

sqrt(2 ^ 3)
sqrt(value1)
sqrt(value2)
help(sqrt)			
?help




# B. Importing Datasets 

getwd()  # Print the current working directory
setwd("C:/Users/Sonke.Ehret/Document/ExperimentCourse")  # Set up a new working directory, 
dir()

# 1. Read .csv format (comma separated values) 
health <- read.csv("Dataset.csv")    
health <- read.table("Dataset.csv", sep = ",", header = TRUE)  
class(health)
# Or equivalently you can include the path--  
# read.csv("/Users/Sonke.Ehret/Document/Experimentcourse")   


# 2. Read .txt format 
# read.table("Dataset.txt", sep = " ")

# C. Packages
library() 			# Check all packages installed
search() 			  # Check pagkages currently loaded

# To Install a New Package
# - Step1: install.packages("name of package")
# - Step2: library(name of package) 

# 2.a To read in text data faster, use the "readr" package
# Save dataset.csv as text file Dataset.txt and use "readr" to re-import
read_table('Dataset.txt') # when your data are separated by one or more spaces
read_delim('Dataset.txt', delim = '\t')
read_csv('Dataset.csv')  


# 3. Read SAS, SPSS, and Stata Format: package "haven"
# install.packages("haven")
# library(haven)
# read_spss("Dataset.sav")
# read_dta("Dataset.dta")
# read_sas("Dataset.sas7bdat)


# 4. Read Excel Format 
# install.packages("readxl")
# package reads both xls and xlsx files
# library(readxl)
# read_excel("Dataset.xlsx")


# D. Export data file 
write.csv(health, "IntrotoR.final.csv", row.names = FALSE) 

# Use the readr package
# write_csv(health, 'healthExam.csv')

# Use the haven package to export SPSS or Stata files
# write_spss(health, "my_spss.sav")
# write_dta(health, "my_stata.dta")

## 

# ----------------------------------------------------------
# -- II. Randomization Inference                        ----
# ----------------------------------------------------------  

## simulation parameters
n <- 100 # sample size
mu0 <- 0 # mean of Y_i(0)
sd0 <- 1 # standard deviation of Y_i(0)
mu1 <- 1 # mean of Y_i(1)
sd1 <- 1 # standard deviation of Y_i(1)

## generate a sample
Y0 <- rnorm(n, mean = mu0, sd = sd0)
Y1 <- rnorm(n, mean = mu1, sd = sd1)
tau <- Y1 - Y0 # individual treatment effect
## true value of the sample average treatment effect
SATE <- mean(tau)
SATE

## repeatedly conduct randomized controlled trials
sims <- 5000 # repeat 5,000 times, we could do more
diff.means <- rep(NA, sims)  # container

for (i in 1:sims) {
  ## randomize the treatment by sampling of a vector of 0's and 1's
  treat <- sample(c(rep(1, n / 2), rep(0, n / 2)), size = n, replace = FALSE)
  ## difference-in-means
  diff.means[i] <- mean(Y1[treat == 1]) - mean(Y0[treat == 0])
}

## estimation error for SATE
est.error <- diff.means - SATE
summary(est.error)



# ----------------------------------------------------------
# -- III. Randomization Inference Using the ri package  ----
# ----------------------------------------------------------  

# simulation in which 2 of 7 villages from Table '2.1' (lecture) are assigned to treatment

rm(list=ls())       # clear objects in memory
library(ri)         # load the RI package
set.seed(1234567)   # random number seed, so that results are reproducible

# input full schedule of potential outcomes
# using Table 2.1

Y0 <- c(10,15,20,20,10,15,15)
Y1 <- c(15,15,30,15,20,15,30)

# create a potential outcomes object called a data frame

Ys <- data.frame(Y0,Y1)
# check column means
colMeans(Ys)

# create a vector of possible assignments 
Z  <- c(rep(1,2),rep(0,5))

# in order to randomly sample with replacement from Z
# type the command
# sample(Z)

# generate all permutations of Z under _complete_ random assignment
# note that default is to do every possible permutation if less than 10,000 permutations

perms <- genperms(Z)

# show number of permutations
cat(ncol(perms)," = number of permutations") 

probs <- genprobexact(Z,blockvar=NULL)  # inputs imply equal-probability assignment
# verify that probability of treatment is constant across the sample
table(probs)

# calculate the sampling distribution of estimated difference-in-means
truedist <- gendist(Ys,perms,Ypre=NULL,prob=probs,HT=FALSE)

# display the frequency distribution of the sampling distribution
table(truedist)

# graphically display the sampling distribution
dispdist(truedist,0)

# show the ATE estimate for each random assignment
truedist



sum((truedist-mean(truedist))^2)

se<-sqrt(sum((truedist-mean(truedist))^2)/length(truedist))


# calculate the proportion of estimates that are above zero

length(truedist[truedist > 0])
length(truedist[truedist > 0])/length(truedist)

# -------------------------------------------
# -- IV. Block Randomization   --
# -------------------------------------------  

rm(list=ls(all=TRUE))
library(ri)
set.seed(1234567)

Y0 <- c(0,1,2,4,4,6,6,9,14,15,16,16,17,18)
Y1 <- c(0,0,1,2,0,0,2,3,12,9,8,15,5,17)

Z <- c(1,1,0,0,0,0,0,0,0,0,0,0,1,1)

# generate all permutations of Z under _complete_ random assignment
# note that default is to do every possible permutation if less than 10,000 permutations

compperms <- genperms(Z)
numperms <- ncol(compperms)

# create empty vector
compmeans <- rep(NA,numperms)

# loop to create average treatment effect estimates for each randomization
for (i in 1:numperms) compmeans[i] <- mean(Y1[compperms[,i]==1]) - mean(Y0[compperms[,i]==0])

# randomize within blocks
block <- c(1,1,1,1,1,1,1,1,2,2,2,2,2,2)

# generate all permutations of Z under block random assignment

blockperms <- genperms(Z,block)
numperms <- ncol(blockperms)

# create empty vector
blockmeans <- rep(NA,numperms)

# loop to create average treatment effect estimates for each randomization
for (i in 1:numperms) blockmeans[i] <- weighted.mean(Y1[blockperms[,i]==1],c(8/2,8/2,6/2,6/2)) - weighted.mean(Y0[blockperms[,i]==0],c(8/6,8/6,8/6,8/6,8/6,8/6,6/4,6/4,6/4,6/4))

save(compmeans,blockmeans,file="figure3.1.Rdata")

# Draw histograms for Figure 3.1	

par(mfrow=c(2,1))
hist(compmeans,main="Sampling Distribution under Complete Randomization",xlim=c(-15,10),xlab="ATE Estimates",freq=FALSE,ylim=c(0,.3))
hist(blockmeans,main="Sampling Distribution under Blocked Randomization",xlim=c(-15,10),xlab="ATE Estimates",freq=FALSE,ylim=c(0,.3))

# calculate the proportion of estimates that are above zero

length(compmeans[compmeans > 0])
length(compmeans[compmeans > 0])/length(compmeans)

length(blockmeans[blockmeans > 0])
length(blockmeans[blockmeans > 0])/length(blockmeans)


# -------------------------------------------
# -- V. Homework   --
# -------------------------------------------  


# 1. Gerber & Green 3.4

# 2. Gerber & Green 3.6
# dataset:
dataset36<-read.csv('http://hdl.handle.net/10079/ghx3frr')  

