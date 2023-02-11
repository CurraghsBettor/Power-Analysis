rm(list=ls()) # clear workspace

# for reproducibility
set.seed(123)

# number of participants per group
n <- 30

library(MASS)
# set the different population means:
mu1 <- c(0.36,0.36)
mu2 <- c(1.4, 1.4)
mu3 <- c(-0.8, 0.8)

# create the variance covariance matrix with Sigma = 1 and Rho = 0.7
cov_matrix <- array(c(1,0.7,0.7,1), dim =c(2,2))
cov_matrix

# generate random data from bivariate normal distributions
gp1 <- mvrnorm(n, mu = mu1, Sigma = cov_matrix)
gp2 <- mvrnorm(n, mu = mu2, Sigma = cov_matrix)
gp3 <- mvrnorm(n, mu = mu3, Sigma = cov_matrix)

# displays the correlation between repeated measures in each group 
cor1 <- cor(x = gp1[,1], y = gp1[,2]); print(cor1)
cor2 <- cor(x = gp2[,1], y = gp2[,2]); print(cor2)
cor3 <- cor(x = gp3[,1], y = gp3[,2]); print(cor3)

# create the data frame
factor1 <- c(gp1[,1], gp2[,1], gp3[,1])
factor2 <- c(gp1[,2], gp2[,2], gp3[,2])

id <- 1:90
group <- c(rep(1, 30), rep(2, 30), rep(3, 30))
data <- cbind(id, group, factor1, factor2)
df <- as.data.frame(data); print(df)

# reshape to long format
library(reshape2)
data_long <- melt(df, id.vars = c("id", "group"), measure.vars = c("factor1", "factor2"), variable.name = "factor", value.name = "DV")

# run the anovas
options(contrasts = c("contr.sum", "contr.poly"))
library(ez)
data_long$id <- as.factor(data_long$id)
data_long$group <- as.factor(data_long$group)
data_long$factor <- as.factor(data_long$factor)
data_long$DV <- as.numeric(data_long$DV)
Model <- ezANOVA(data_long,
                 dv =. (DV), wid =.(id),
                 within =. (factor),
                 between =. (group),
                 type = 3, 
                 detailed = T,
                 return_aov = T)
report(Model)

Model2 <- ezANOVA(data_long,
                  dv =. (DV), wid =.(id),
                  between =. (group),
                  type = 3,
                  detailed = T,
                  return_aov = T); report(Model2)

# estimate the correlation between repeated measures regardless of the group
corG <- cor(x = factor1, y = factor2); print(corG)

Model3 <- ezANOVA(data_long,
                  dv =. (DV), wid =.(id),
                  within =. (factor),
                  type = 3,
                  detailed = T,
                  return_aov = T); report(Model3)

# estimate partial etas when factor variable is considered as both within and between-subjects variables
# compute the t-statistics
tb <- t.test(x = factor1, y = factor2, paired = F, var.equal = T)$statistic
tw <- t.test(x = factor1, y = factor2, paired = T)$statistic
# find each corresponding partial eta squared
etab <- tb^2/sum(tb^2, 178); print(etab)
etaw <- tw^2/sum(tw^2, 89); print(etaw)

# set the parameters n = 90, Delta = 15, Sigma = 15, and Rho = 0.5 and Rho = 0
# Compute the different Cohen's d

Delta <- 15; Sigma <- 15; Rho <- c(0, 0.5); n = 90

sdp <- sqrt((Sigma^2*(n-1)+Sigma^2*(n-1))/(sum(n,n)-2))

sdav <- sqrt(sum(Sigma^2, Sigma^2)/2)

sdD <- c()
for(i in Rho) {
  stand <- sqrt(Sigma^2+Sigma^2-2*Sigma*Sigma*i)
  sdD <- c(sdD, stand)
}

# Cohen’s d_p
dp <- Delta/sdp; dp

# Cohen’s d_av = Cohen’s d_p
dav <- Delta/sdav; dav

# Cohen’s d_D with Rho = 0.5
dDHalf <- Delta/sdD[2]; dDHalf

# Cohen’s d_D with Rho = 0
dDZero <- Delta/sdD[1]; dDZero

library(psych)
ms <- c(n, n)
hm <- harmonic.mean(ms)

## compute the t-statistic/lambda
# between-subjects design
lambda_dp <- dp* sqrt(hm/2)
# within-subjects design, when Rho = 0.5
lambda_dDH <- dDHalf*sqrt(n)
# within-subjects design, when Rho = 0
lambda_dDZ <- dDZero*sqrt(n)
## compute the degrees of freedom
dfp <- sum(n,n)-2
dfD <- n-1

# function to obtain a partial eta squared from the t-statistic and the degrees of freedom
eta <- function(t, df) {
  partial <- t^2/sum(t^2, df)
  return(partial)
}

# eta for a between subjects-variable that is equivalent to the generalized eta squared
etaB <- eta(lambda_dp, dfp); eta
# eta for a within subjects-variable when Rho = 0.5
etawH <- eta(lambda_dDH, dfD); etawH
# eta for a within subjects-variable when Rho = 0
etawZ <- eta(lambda_dDZ, dfD); etawZ

# function that allows to convert a partial eta squared for a within-subjects variable as in SPSS to a partial eta squared as in GPower (i.e.,that is comparable to a partial eta squared computed for a between-subjects variable) 

Spss_to_GP <-function(eta, N, k, m, Rho) { 
  fsq <- eta/(1-eta)
  Gpower <-fsq*((N-k)/N)*((m-1)/m)*(1-Rho)
  etaGP <-GPower/sum(1, 1, GPower)
  return(test)
}

