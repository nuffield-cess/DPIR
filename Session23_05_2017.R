
# Remove every objects in your working environment 
rm(list = ls())


# -------------------------------------------
# -- I.Mediation   --
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

med.out <- ?mediate(med.fit, out.fit, treat = "treat", mediator = "emo",
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







