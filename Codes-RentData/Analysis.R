
library(Rfast)
library(hdm)
library(glmnet)
library(sandwich)
library(lmtest)
library(stargazer)
library(xtable)
library(ggplot2)

source("Functions.R")

# FILE WHERE ALL METHODS START THERE CALCULATIONS

Final <- read.csv("/Users/franziskaaupperle/Desktop/Finals/Rent-Data/Data-full/Final.csv")
Final_dmy <- read.csv("/Users/franziskaaupperle/Desktop/Finals/Rent-Data/Data-Full/Final_dmy.csv")
Final_dmy1 <- read.csv("/Users/franziskaaupperle/Desktop/Finals/Rent-Data/Data-Full/Final_dmy1.csv")
Final <- Final[,-1]
Final_dmy <- Final_dmy[,-1]
Final_dmy1 <- Final_dmy1[,-1]

# initialize variables
Rent <- Final_dmy1[,1]
Rental.Brake <- Final_dmy1[,2] 
Controls <- as.matrix(Final_dmy1[,-c(1:2)])

y <- Rent
d <- Rental.Brake
x <- Controls


# PLAIN OLS ESTIMATION
OLS <- lm(y ~ d + x)
# OLS with robust SE
OLS_white <- coeftest(OLS, vcov = vcovHC(OLS, type="HC1"))
# Output TeX Code
stargazer(OLS,
          OLS_white,
          no.space = TRUE, save = TRUE, single.row = TRUE)

# prepare data for inference on rental brake with proposed methods

# DOUBLE-SELECTION
Double <- rlassoEffect(y = y, d = d, x = x, method = "double selection", I3 = NULL)
# get selection indices
DB <- Double$selection.index
# refit the model
lm.Double <- lm(y ~ d + x[,DB])
# refit the model with robust standard errors
lm.Double_white <- coeftest(lm.Double, vcov = vcovHC(lm.Double, type="HC1"))
stargazer(lm.Double,
          lm.Double_white,
          no.space = TRUE, save = TRUE, single.row = TRUE)



# PODS
PODS <- Projection(y = y, d = d, x = x)
# get selection indices
PDS <- PODS$M.hat
# refit the model
lm.PODS <- lm(y ~ d + x[,PDS])
# refit the model with robust standard errors
lm.PODS_white <- coeftest(lm.PODS, vcov = vcovHC(lm.PODS, type="HC1"))
stargazer(lm.PODS,
          lm.PODS_white,
          no.space = TRUE, save = TRUE, single.row = TRUE)


# note: whether LASSO in the upcoming 3 methods is performed with cross-validation or the rigorous LASSO, 
# needs to be defined in function select.model()

# R-SPLIT
# perform R-Split with LASSO model selection applying cross-validation
RSPLIT_W_CV <- Split.smooth(y = y, d = d, x = x, B = 1000)
# perform R-Split with LASSO model selection applying rigorous LASSO
RSPLIT_WO_CV <- Split.smooth(y = Rent, d = d, x = x, B = 1000)


# PODS-SPLIT
# perform PODS-Split with LASSO model selection applying cross-validation
PODS.SPLIT_W_CV <- PODS.Split(y = y, d = d, x = x, B = 1000)
# perform PODS-Split with LASSO model selection applying rigorous LASSO
PODS.SPLIT_WO_CV <- PODS.Split(y = y, d = d, x = x, B = 1000)


# DOUBLE STABILITY
# perform Double-Stability with LASSO model selection applying cross-validation
DB.STABILITY_W_CV <- Double.Stability(y = y, d = d, x = x, B = 1000)
# perform Double-Stability with LASSO model selection applying rigorous LASSO
DB.STABILITY_WO_CV <- Double.Stability(y = y, d = d, x = x, B = 1000)



# Prepare data for analysis of relative incl frequencies
# Analysis can only run for one method of R-Split, PODS-SPlit and Double-Stability

# R-Split - Select relative inclusion frequencies for each variable from either the cross-validated approach or the rigorous Lasso
# freq <- colMeans(RSPLIT_W_CV[["C.count"]])
# freq <- colMeans(RSPLIT_WO_CV[["C.count"]])

# PODS.Split - Select relative inclusion frequencies for each variable from either the cross-validated approach or the rigorous Lasso
# freq <- colMeans(PODS.SPLIT_W_CV[["C.count"]])
# freq <- colMeans(PODS.SPLIT_W0_CV[["C.count"]])

# Double-Stability - Select relative inclusion frequencies for each variable from either the cross-validated approach or the rigorous Lasso
df <- rbind(colMeans(DB.STABILITY_W_CV[["Cd.count"]]),  colMeans(DB.STABILITY_W_CV[["Cy.count"]]))
# df <- rbind(colMeans(DB.STABILITY_WO_CV[["Cd.count"]]),  colMeans(DB.STABILITY_WO_CV[["Cy.count"]]))
freq <- round(apply(df, 2, max),2)

# declare data frame for each variable which depends on the choice of the treshold/relative inclusion frequency
RB.est <- data.frame(data = NA, nrow = length(freq), ncol = 2)
RB.se <- data.frame(data = NA, nrow = length(freq), ncol = 2)
model.size <- data.frame(data = NA, nrow = length(freq), ncol = 2)
R2 <- data.frame(data = NA, nrow = length(freq), ncol = 2)
R2_adj <- data.frame(data = NA, nrow = length(freq), ncol = 2)

for (i in 1:length(freq)) {
  
  # get index of the variable with the i smallest relative inculsion frequency
  ndx <- order(freq)[i]
  print(ndx)
  print(freq[ndx])
  # get all inidices that have a relative inclusion frequency of at least freq[ndx]
  M.hat <- which(freq >= freq[ndx])
  print(colnames(x[,M.hat]))
  
  # Post-LASSO OLS step
  # get estimator for the rental brake
  lm <- lm(y ~ d + x[,M.hat]) 
  RB.est[i,1] <- freq[ndx]
  RB.est[i,2] <- lm$coef[2]
  
  # get its HC1 standard error
  lm_HC1 = coeftest(lm, vcov = vcovHC(lm, type="HC1"))
  RB.se[i,1] <- freq[ndx]
  RB.se[i,2] <- lm_HC1[2,2]
 
  # get its model size
  model.size[i,1] <- freq[ndx]
  model.size[i,2] <- length(M.hat)

  # get R-Squared of the regression
  R2[i,1] <- freq[ndx]
  R2[i,2] <- summary(lm)$r.squared

  # get its Adj. R2        
  R2_adj[i,1] <- freq[ndx]
  R2_adj[i,2] <- summary(lm)$adj.r.squared
  
}

# Grafics for the Estimated effects, Adj. R^2 and Model Size depending on the treshold
# As well as Empirical Cumulative Distribution Function for the realtive inclusion frequ in R-Split and PODS-Split 

# ECDF for the RIF
ggplot(model.size, aes(model.size[,1])) + stat_ecdf(geom = "step")+
  labs(title="ECDF for RIF \n of Variables in R-Split",
       y = "Cumulative Frequency", x="Relative Inclusion Frequency") +
  theme_classic()

# Plot for the Rental Brake depending on the treshold
ggplot() + geom_line(RB.est, mapping = aes(x = RB.est[,1], y = RB.est[,2]))+
  labs(title="Estimated Coefficient of Rental.Brake \n in the Refitted Model",
       y = "Estimated coefficient of 'Rental.Brake'", x="Treshold") +
  theme_classic()

# Plot for the Adj. R2 depending on the treshold
ggplot() + geom_line(R2_adj, mapping = aes(x = R2_adj[,1], y = R2_adj[,2]))+
  labs(title="Adjusted R-Squared \n of the Refitted Model",
       y = "Adjusted R-Squared", x="Treshold") +
  theme_classic()

# Plot for the model size depending on the treshold
ggplot() + geom_line(model.size, mapping = aes(x = model.size[,1], y = model.size[,2]))+
  labs(title="Model Size of the Refitted Model",
       y = "Model Size", x="Treshold") +
  theme_classic()










