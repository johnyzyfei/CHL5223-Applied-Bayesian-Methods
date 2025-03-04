# Load Required Libraries
library(rjags)
library(coda)  # For MCMC diagnostics and summaries
# Author: Dr. Michael Escobar 


# Description of the dataset
cat("
/*   
Data taken from page 204 of McCullagh Nelder (2nd ed)

The data was provided by J. Drilley and L.N. Hemingway of Lloyd's
Register of Shipping, concern a type of damage caused by waves to the forward 
section of certain cargo-carrying vessels. For the purpose of setting standards for hull 
construction we need to know the risk of damage associated with the three
classifying factors shown below:

Column 1:  Ship type, 1-5
Column 2: Year of construction: 1=1960-64, 2=1965-69, 3=1970-74, 4=1975-79
Column 3: Period of operation: 1=1960-74, 2=1975-1979
Column 4: Aggregate months of service
Column 5: Number of reported damage incidents
*/", file="WaveDamDesc.txt")

# Write Data to File
cat(
  "ship   yrcons   yrop  month daminc
1   1  1   127     0
1   1  2   63      0
1   2  1   1095    3
1   2  2   1095    4
1   3  1   1512    6
1   3  2   3353    18
1   4  2   2244    11
2   1  1  44882    39
2   1  2  17176    29
2   2  1  28609    58
2   2  2  20370    53
2   3  1   7064    12
2   3  2  13099    44
2   4  2   7117    18
3   1  1   1179    1
3   1  2    552    1
3   2  1    781    0
3   2  2    676    1
3   3  1    783    6
3   3  2   1948    2
3   4  2    274    1
4   1  1    251    0
4   1  2    105    0
4   2  1    288    0
4   2  2    192    0
4   3  1    349    2
4   3  2   1208    11
4   4  2   2051    4
5   1  1     45    0
5   2  1    789    7
5   2  2    437    7
5   3  1   1157    5
5   3  2   2161    12
5   4  2    542    1
", file="WaveDamDat.txt")

# Define Bayesian Model for JAGS
cat("
model {
  for(i in 1:34) {
    daminc[i] ~ dpois(lam[i])
    log(lam[i]) <- log(month[i]) + beta0 + beta.o*yrop[i] + beta.s[ship[i]] + beta.c[yrcons[i]] + b[i] 
    b[i] ~ dnorm(0, tau)
    b.adj[i] <- b[i] - mean(b[])
  }
  
  # Ship type effects
  for(is in 1:5) {
    beta.s[is] ~ dnorm(0, tau.s)
    beta.s.adj[is] <- beta.s[is] - mean(beta.s[])
  }
  
  # Construction year effects
  for(ic in 1:4) {
    beta.c[ic] ~ dnorm(0, tau.c)
    beta.c.adj[ic] <- beta.c[ic] - mean(beta.c[])
  }

  # Prior for base rate
  beta0 ~ dnorm(0, .0082)
  beta0.adj <- beta0 + mean(b[]) + mean(beta.s[]) + mean(beta.c[])

  # Extra Poisson variation
  std ~ dunif(0, 9)
  tau <- 1 / std / std

  # Priors for relative risk
  std.s ~ dunif(0, 5)
  tau.s <- 1 / std.s / std.s
  std.c ~ dunif(0,5)
  tau.c <- 1 / std.c / std.c

  beta.o ~ dnorm(0, .04)
}
", file="waveModCentered.txt")

# Load Data
WaveDamDat <- read.table("WaveDamDat.txt", header=TRUE, sep = "")
bug.dat <- list(
  ship = WaveDamDat$ship,
  yrcons = WaveDamDat$yrcons,
  yrop = WaveDamDat$yrop,
  month = WaveDamDat$month,
  daminc = WaveDamDat$daminc
)

# Initial Values Function
init.fun <- function() list(
  beta.s = rnorm(5), std.s = runif(1,1,2),
  beta.c = rnorm(4), std.c = runif(1,1,2), 
  beta.o = rnorm(1), std = runif(1,1,2), beta0 = rnorm(1),
  b = rnorm(34,0,.1)
)

# Parameters to Monitor
params <- c("beta.s.adj", "std.s", "beta.c.adj", "std.c", "beta.o", "beta0.adj", "std")

# Run JAGS Model
jags_model <- jags.model("waveModCentered.txt", data = bug.dat, inits = init.fun, n.chains = 5, n.adapt = 500)
update(jags_model, 1000)  # Burn-in

# Sample from posterior
WaveBug0 <- coda.samples(jags_model, variable.names = params, n.iter = 2000)

# Print Summary Results
summary(WaveBug0)

# Boxplots for Ship and Construction Effects
par(mfrow=c(1,2))
boxplot(as.matrix(WaveBug0[[1]])[,grep("beta.s.adj", colnames(as.matrix(WaveBug0[[1]])))], main="Ship Type Effects (beta.s.adj)")
boxplot(as.matrix(WaveBug0[[1]])[,grep("beta.c.adj", colnames(as.matrix(WaveBug0[[1]])))], main="Construction Year Effects (beta.c.adj)")

# Autocorrelation Analysis
par(mfrow=c(2,2))
acf(as.matrix(WaveBug0[[1]])[, "beta.s.adj[1]"], main="ACF: beta.s.adj[1]")
acf(as.matrix(WaveBug0[[1]])[, "beta.c.adj[1]"], main="ACF: beta.c.adj[1]")
acf(as.matrix(WaveBug0[[1]])[, "beta0.adj"], main="ACF: beta0.adj")
acf(as.matrix(WaveBug0[[1]])[, "beta.o"], main="ACF: beta.o")

# Density Plots
par(mfrow=c(2,2))
plot(density(as.matrix(WaveBug0[[1]])[, "beta.s.adj[2]"]), main="Density: beta.s.adj[2]")
plot(density(as.matrix(WaveBug0[[1]])[, "beta.c.adj[2]"]), main="Density: beta.c.adj[2]")

# Trace Plots
matplot(as.matrix(WaveBug0[[1]]), type="l", main="Trace Plots")

