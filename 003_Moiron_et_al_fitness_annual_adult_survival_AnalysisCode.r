# Code for "The quantitative genetics of fitness in a wild bird population"
# Unpublished manuscript, doi: tba
# Moiron M, Charmantier A, Bouwhuis S

# The code provided here is sufficient to replicate the simulations presented in the above paper

######################################################
# DATA ANALYSIS OF ANNUAL ADULT SURVIVAL DATA
######################################################

# Loading packages
library(tidyverse)
library(MCMCglmm)
library(QGglmm)
library(pedantics)
library(beepr)

# Loading phenotypic data
data <- read.table("fitness_data.txt", header=TRUE)

# Formating response variable
data$AAS<-as.numeric(data$AAS)  #annual adult survival data

# Formating random effects
data$ID.PE =as.factor(data$ring.nr) #permanent environmental effects
data$animal = as.factor(data$ring.nr) #additive genetic effects
data$YEAR = as.factor(data$year) #year effects

# Formating fixed effects
data$age=as.numeric(scale(data$age)) #linear component of age

# Loading pedigree
pedigree <- read.table("pedigree.txt",header=TRUE)
my_inverse <- inverseA(pedigree)$Ainv

# Setting number of samples and iterations
nsamp <- 1000 
BURN <- 10000; THIN <- 10000
(NITT <- BURN + THIN*nsamp)

# Setting prior
priorAAS <- list(R = list(V = 1, fix = 1),
              G = list(G1 = list(V = 1, nu = 2, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(100)),
                       G2 = list(V = 1, nu = 2, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(100)),
                       G3 = list(V = 1, nu = 2, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(100))))

# Running model
mod<- MCMCglmm(AAS ~ 1 + age,
               random = ~ animal +ID.PE +YEAR,
               data = data, 
               ginverse = list(animal = my_inverse),
               prior = priorAAS, 
               family = "categorical",
               nitt = NITT, thin = THIN, burnin = BURN, 
               verbose = TRUE)

beep(1)
#save(mod, file = "analysis_ARS.rda")
#load("analysis_ARS.rda")

summary(mod)

# Assesing convergence and auto-correlation of posterior distributions
plot(mod$VCV)
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)
autocorr.diag(mod$VCV)

# Random effects in latent scale
round(posterior.mode(mod$VCV),3)
round(HPDinterval(mod$VCV), 3)

# Computing h2 and Ia
df <- data.frame(mu = as.vector(mod[["Sol"]][ , "(Intercept)"]),
                 va = as.vector(mod[["VCV"]][ , "animal"]), vp = rowSums(mod[["VCV"]]))
post <-  do.call("rbind",apply(df, 1, function(row) {
  QGparams(mu = row["mu"], var.a = row["va"], var.p = row["vp"],
           model = "binom1.probit", verbose = FALSE)
}))

post[["Ia.obs"]] <- post[["var.a.obs"]] / (post[["mean.obs"]]^2)

# Saving output
out<- data.frame(mean.obs.mode = posterior.mode(as.mcmc(post[["mean.obs"]])),
                 mean.obs.mean = mean(as.mcmc(post[["mean.obs"]])),
                 mean.obs.low  = HPDinterval(as.mcmc(post[["mean.obs"]]))[1],
                 mean.obs.high = HPDinterval(as.mcmc(post[["mean.obs"]]))[2],
                 va.mode       = posterior.mode(as.mcmc(post[["var.a.obs"]])),
                 va.mean       = mean(as.mcmc(post[["var.a.obs"]])),
                 va.low        = HPDinterval(as.mcmc(post[["var.a.obs"]]))[1],
                 va.high       = HPDinterval(as.mcmc(post[["var.a.obs"]]))[2],
                 vp.mode       = posterior.mode(as.mcmc(post[["var.obs"]])),
                 vp.mean       = mean(as.mcmc(post[["var.obs"]])),
                 vp.low        = HPDinterval(as.mcmc(post[["var.obs"]]))[1],
                 vp.high       = HPDinterval(as.mcmc(post[["var.obs"]]))[2],
                 h2.mode       = posterior.mode(as.mcmc(post[["h2.obs"]])),
                 h2.mean       = mean(as.mcmc(post[["h2.obs"]])),
                 h2.low        = HPDinterval(as.mcmc(post[["h2.obs"]]))[1],
                 h2.high       = HPDinterval(as.mcmc(post[["h2.obs"]]))[2],
                 Ia.mode       = posterior.mode(as.mcmc(post[["Ia.obs"]])),
                 Ia.mean       = mean(as.mcmc(post[["Ia.obs"]])),
                 Ia.low        = HPDinterval(as.mcmc(post[["Ia.obs"]]))[1],
                 Ia.high       = HPDinterval(as.mcmc(post[["Ia.obs"]]))[2])
out
