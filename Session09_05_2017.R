library(randomizr)


# -------------------------------------------
# -- I. Cluster random assignment  --
# -------------------------------------------  

#A. Setup

# Load built-in dataset
data(HairEyeColor)
HairEyeColor <- data.frame(HairEyeColor)

# Transform so each row is a subject
# Columns describe subject's hair color, eye color, and gender
hec <- HairEyeColor[rep(1:nrow(HairEyeColor),
                        times = HairEyeColor$Freq), 1:3]

N <- nrow(hec)


# Fix the rownames
rownames(hec) <- NULL

# Set a seed for reproducability
set.seed(343)

# Create untreated and treated outcomes for all subjects
hec <- within(hec,{
  Y0 <- rnorm(n = N,mean = (2*as.numeric(Hair) + -4*as.numeric(Eye) + -6*as.numeric(Sex)), sd = 5)
  Y1 <- Y0 + 6*as.numeric(Hair) + 4*as.numeric(Eye) + 2*as.numeric(Sex)
})

# Calculate true ATE
with(hec, mean(Y1 - Y0))


# B. Simple Random assignment

Z <- simple_ra(N = N)
table(Z)
Z <- simple_ra(N = N, prob = 0.30)

Z <- simple_ra(N = N, num_arms = 3)

Z <- simple_ra(N = N, prob_each = c(.2, .2, .6))


Z <- simple_ra(N = N, prob_each = c(.2, .2, .6),
               condition_names=c("control", "placebo", "treatment"))
table(Z)

#C. Complete random assignment

Z <- complete_ra(N = N)
table(Z)

Z <- ?complete_ra(N = N, m = 200)
table(Z)

Z <- complete_ra(N = N, num_arms = 3)
table(Z)

Z <- complete_ra(N = N, m_each = c(100, 200, 292),
                 condition_names = c("control", "placebo", "treatment"))
table(Z)

# C. Block random assignment 

Z <- block_ra(block_var = hec$Hair)
table(Z, hec$Hair)


Z <- block_ra(block_var = hec$Hair, num_arms = 3)
table(Z, hec$Hair)


Z <- block_ra(block_var = hec$Hair, condition_names = c("Control", "Placebo", "Treatment"))
table(Z, hec$Hair)


Z <- block_ra(block_var = hec$Hair, prob_each = c(.3, .7))
table(Z, hec$Hair)


sort(unique(hec$Hair))


block_m_each <- rbind(c(78, 30),
                      c(186, 100),
                      c(51, 20),
                      c(87,40))


block_m_each


Z <- block_ra(block_var = hec$Hair, block_m_each = block_m_each)
table(Z, hec$Hair)

#Try with


block_m_each_alt <- rbind(c(78, 20, 10),
                          c(186, 90, 10),
                          c(51, 10, 10),
                          c(87,30, 10))


declaration <- ?declare_ra(block_var = hec$Hair, block_m_each = block_m_each)
cond_prob <- obtain_condition_probabilities(declaration, Z)
table(cond_prob, Z)
# show the probability that each unit is assigned to each condition
head(declaration$probabilities_matrix)


# Show that the probability of treatment is different within block
table(hec$Hair, round(declaration$probabilities_matrix[,2], 3))


hec <- within(hec,{
  Z_blocked <- block_ra(block_var = hec$Hair,
                        block_m_each = block_m_each)
  Y_blocked <- Y1*(Z_blocked) + Y0*(1-Z_blocked)
  cond_prob <- obtain_condition_probabilities(declaration, Z_blocked)
  IPW_weights <- 1/(cond_prob)
})

fit_LSDV <- lm(Y_blocked ~ Z_blocked + Hair, data=hec)
fit_IPW <- lm(Y_blocked ~ Z_blocked, weights = IPW_weights, data = hec)

summary(fit_LSDV)

summary(fit_IPW)

block_var <- with(hec, paste(Hair, Eye, Sex, sep = "_"))
Z <- block_ra(block_var = block_var)
head(table(block_var, Z))


library(blockTools)

# BlockTools requires that all variables be numeric
numeric_mat <- model.matrix(~Hair+Eye+Sex, data=hec)[,-1]

# BlockTools also requres an id variable
df_forBT <- data.frame(id_var = 1:nrow(numeric_mat), numeric_mat)

# Conducting the actual blocking: let's make trios
out <- block(df_forBT, n.tr = 3, id.vars = "id_var", 
             block.vars = colnames(df_forBT)[-1])

# Extact the block_ids
hec$block_id <- createBlockIDs(out, df_forBT, id.var = "id_var")

# Conduct actual random assignment with randomizr
Z_blocked <- block_ra(block_var = hec$block_id, num_arms = 3)
head(table(hec$block_id, Z_blocked))


# D. Clustered Assignment

clust_var <- with(hec, paste(Hair, Eye, Sex, sep = "_"))
hec$clust_var <- clust_var

Z_clust <- cluster_ra(clust_var = clust_var)

head(table(clust_var, Z_clust))

Z_clust <- cluster_ra(clust_var = clust_var, num_arms = 3)
head(table(clust_var, Z_clust))

Z_clust <- cluster_ra(clust_var=clust_var, 
                      condition_names=c("Control", "Placebo", "Treatment"))
head(table(clust_var, Z_clust))


Z_clust <- cluster_ra(clust_var=clust_var, m_each=c(5, 15, 12))
head(table(clust_var, Z_clust))



# -------------------------------------------
# -- II. Regressions & Covariates   --
# ------------------------------------------- 


#Clear any previous work
rm(list=ls(all=TRUE))

#Load Relevant packages
library(AER)
library(sandwich)

data1 <- read.csv(file="http://hdl.handle.net/10079/70rxwqn",head=TRUE,sep=",")

# You can check that your data was read in correctly using these two commands:
# colnames(data1)
# dim(data1)

# select one-person households that were either pure controls or canvass only
sel <-  data1$onetreat==1 & data1$mailings==0 & data1$phongotv==0 & data1$persons==1

# verify the number of observations
table(sel)
data2 <- data1[sel,]

# define variables
v98      <- data2$v98
persngrp <- data2$persngrp
cntany   <- data2$cntany

############  NOTE USE OF ROBUST STANDARD ERRORS


# Box 5.4: ITT
coef(summary(lm(v98 ~ persngrp)))
# robust SEs
itt_fit <- lm(v98 ~ persngrp)
coeftest(itt_fit,vcovHC(itt_fit))


# Box 5.5: ITT_D
# Note that results from this will vary based on the current version that you have but this variation should not be a concern. 
coef(summary(lm(cntany ~ persngrp)))
# robust SEs
itt_d_fit <- lm(cntany ~ persngrp)
coeftest(itt_d_fit,vcovHC(itt_d_fit))


# Box 5.6: CACE
coef(summary(ivreg(v98 ~ cntany,~persngrp)))
# robust SEs
cace_fit <- ivreg(v98 ~ cntany,~persngrp)
coeftest(cace_fit,vcovHC(cace_fit))


# Program to create the regression boxes for Chapter 6 using NYC Debates data  Mullainathan et al 2010

#Clear any previous work
rm(list=ls(all=TRUE))

library(foreign)
library(AER)



# -------------------------------------------
# -- III.Mediation   --
# ------------------------------------------- 

# Download data online. This is a simulated dataset for this post.
myData <- read.csv('http://static.lib.virginia.edu/statlab/materials/data/mediationData.csv')


model.0 <- lm(Y ~ X, myData)
summary(model.0)

model.M <- lm(M ~ X, myData)
summary(model.M)

model.Y <- lm(Y ~ X + M, myData)
summary(model.Y)

library(mediation)
results <- mediate(model.M, model.Y, treat='X', mediator='M',
                   boot=TRUE, sims=500)
summary(results)

model.M <- lm(M ~ X, myData)
model.Y <- lm(Y ~ X + M, myData)
results <- mediate(model.M, model.Y, treat='X', mediator='M',
                   boot=TRUE, sims=100)
summary(results)

library("mediation")
set.seed(2014)
data("framing", package = "mediation")

?framing

# We use the linear regression fit with least squares and the probit regression
# for the mediator and outcome models, respectively.


med.fit <- lm(emo ~ treat + age + educ + gender + income, data = framing)
out.fit <- glm(cong_mesg ~ emo + treat + age + educ + gender + income,
                  data = framing, family = binomial("probit"))

# We now use the mediate function to estimate the average causal mediation effects ACME and 
# average direct effects ADE.

med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo",
                      robustSE = TRUE, sims = 100)
summary(med.out)

# Treatment mediator interaction. The ACME may take
# different values depending on the baseline treatment status

med.out <- mediate(med.fit, out.fit, boot = TRUE, treat = "treat",
                      mediator = "emo", sims = 100)
summary(med.out)

med.fit <- lm(emo ~ treat + age + educ + gender + income, data=framing)
out.fit <- glm(cong_mesg ~ emo * treat + age + educ + gender + income,
                  data = framing, family = binomial("probit"))
med.out <- mediate(med.fit, out.fit, treat = "treat", mediator = "emo",
                      robustSE = TRUE, sims = 100)
summary(med.out)

test.TMint(med.out, conf.level = .95)





