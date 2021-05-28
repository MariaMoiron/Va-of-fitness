# Code for "The quantitative genetics of fitness in a wild bird population"
# Unpublished manuscript, doi: tba
# Moiron M, Charmantier A, Bouwhuis S

# The code provided here is sufficient to replicate the analyses presented in the above paper.

######################################################
# DATA ANALYSIS OF LIFETIME FITNESS DATA
######################################################

# Loading packages
library(tidyverse)
library(MCMCglmm)
library(QGglmm)
library(pedantics)
library(beepr)

# Loading phenotypic data
data <- read.table("fitness_data.txt", header=TRUE)
data$fitness<-as.numeric(data$LRS_fl) #lifetime fitness data
data$status = as.factor(data$status)  #two levels: alive or dead

# Loading pedigree
pedigree <- read.table("pedigree.txt",header=TRUE)
my_inverse <- inverseA(pedigree)$Ainv

# Setting number of samples and iterations
nsamp <- 1000 
BURN <- 10000; THIN <- 10000
(NITT <- BURN + THIN*nsamp)

# Setting prior
prior <- list(R = list(V = diag(2)/2, nu = 1, fix = 2),
              G = list(G1 = list(V = diag(2)/2, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                       G2 = list(V = diag(2)/2, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000)))

# Running model
mod<- MCMCglmm(fitness ~ trait-1+status,
               random = ~ idh(trait):animal +
                          idh(trait):year,
               rcov = ~ idh(trait):units,
               data = data, 
               ginverse = list(animal = my_inverse), 
               prior = prior, 
               family = "zipoisson",
               nitt = NITT, thin = THIN, burnin = BURN, 
               verbose = TRUE) 

beep(1)
#save(mod, file = "analysis_lifetime_fitness.rda")
#load("analysis_lifetime_fitness.rda")

summary(mod)

# Assesing convergence and auto-correlation of posterior distributions 
plot(mod$VCV)
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)
autocorr.diag(mod$VCV)

# Random effects in latent scale
round(posterior.mode(mod$VCV),3)
round(HPDinterval(mod$VCV), 3)

# Computing h2 and Ia for ZI component
df_zi <- data.frame(mu = as.vector(mod[["Sol"]][ , "traitzi_fitness"]),
                    va = as.vector(mod[["VCV"]][ , "traitzi_fitness.animal"]),
                    vp = rowSums(mod[["VCV"]][ , grepl("traitzi_fitness", colnames(mod[["VCV"]]))]))

post_zi <-  do.call("rbind",apply(df_zi, 1, function(row) {
  QGparams(mu = row["mu"], var.a = row["va"], var.p = row["vp"],
           model = "binom1.logit", verbose = FALSE)
}))

post_zi[["Ia.obs"]] <- post_zi[["var.a.obs"]] / (post_zi[["mean.obs"]]^2)

# Computing h2 and Ia for Poisson component
df_pois <- data.frame(mu = as.vector(mod[["Sol"]][ , "traitfitness"]),
                      va = as.vector(mod[["VCV"]][ , "traitfitness.animal"]),
                      vp = rowSums(mod[["VCV"]][ , grepl("traitfitness", colnames(mod[["VCV"]]))]))

post_pois <-  do.call("rbind",apply(df_pois, 1, function(row) {
  QGparams(mu = row["mu"], var.a = row["va"], var.p = row["vp"],
           model = "Poisson.log", verbose = FALSE)
}))

post_pois[["Ia.obs"]] <- post_pois[["var.a.obs"]] / (post_pois[["mean.obs"]]^2)

# Saving output from both components
out <- data.frame(mean_pois.mode      = posterior.mode(as.mcmc(post_pois[["mean.obs"]])),
                  mean_pois.mean      = mean(as.mcmc(post_pois[["mean.obs"]])),
                  mean_pois.low       = HPDinterval(as.mcmc(post_pois[["mean.obs"]]))[1],
                  mean_pois.high      = HPDinterval(as.mcmc(post_pois[["mean.obs"]]))[2],
                  va_pois.mode        = posterior.mode(as.mcmc(post_pois[["var.a.obs"]])),
                  va_pois.mean        = mean(as.mcmc(post_pois[["var.a.obs"]])),
                  va_pois.low         = HPDinterval(as.mcmc(post_pois[["var.a.obs"]]))[1],
                  va_pois.high        = HPDinterval(as.mcmc(post_pois[["var.a.obs"]]))[2],
                  vp_pois.mode        = posterior.mode(as.mcmc(post_pois[["var.obs"]])),
                  vp_pois.mean        = mean(as.mcmc(post_pois[["var.obs"]])),
                  vp_pois.low         = HPDinterval(as.mcmc(post_pois[["var.obs"]]))[1],
                  vp_pois.high        = HPDinterval(as.mcmc(post_pois[["var.obs"]]))[2],
                  h2_pois.mode        = posterior.mode(as.mcmc(post_pois[["h2.obs"]])),
                  h2_pois.mean        = mean(as.mcmc(post_pois[["h2.obs"]])),
                  h2_pois.low         = HPDinterval(as.mcmc(post_pois[["h2.obs"]]))[1],
                  h2_pois.high        = HPDinterval(as.mcmc(post_pois[["h2.obs"]]))[2],
                  Ia_pois.mode        = posterior.mode(as.mcmc(post_pois[["Ia.obs"]])),
                  Ia_pois.mean        = mean(as.mcmc(post_pois[["Ia.obs"]])),
                  Ia_pois.low         = HPDinterval(as.mcmc(post_pois[["Ia.obs"]]))[1],
                  Ia_pois.high        = HPDinterval(as.mcmc(post_pois[["Ia.obs"]]))[2],
                  mean_zi.mode        = posterior.mode(as.mcmc(post_zi[["mean.obs"]])),
                  mean_zi.mean        = mean(as.mcmc(post_zi[["mean.obs"]])),
                  mean_zi.low         = HPDinterval(as.mcmc(post_zi[["mean.obs"]]))[1],
                  mean_zi.high        = HPDinterval(as.mcmc(post_zi[["mean.obs"]]))[2],
                  va_zi.mode          = posterior.mode(as.mcmc(post_zi[["var.a.obs"]])),
                  va_zi.mean          = mean(as.mcmc(post_zi[["var.a.obs"]])),
                  va_zi.low           = HPDinterval(as.mcmc(post_zi[["var.a.obs"]]))[1],
                  va_zi.high          = HPDinterval(as.mcmc(post_zi[["var.a.obs"]]))[2],
                  vp_zi.mode          = posterior.mode(as.mcmc(post_zi[["var.obs"]])),
                  vp_zi.mean          = mean(as.mcmc(post_zi[["var.obs"]])),
                  vp_zi.low           = HPDinterval(as.mcmc(post_zi[["var.obs"]]))[1],
                  vp_zi.high          = HPDinterval(as.mcmc(post_zi[["var.obs"]]))[2],
                  h2_zi.mode          = posterior.mode(as.mcmc(post_zi[["h2.obs"]])),
                  h2_zi.mean          = mean(as.mcmc(post_zi[["h2.obs"]])),
                  h2_zi.low           = HPDinterval(as.mcmc(post_zi[["h2.obs"]]))[1],
                  h2_zi.high          = HPDinterval(as.mcmc(post_zi[["h2.obs"]]))[2],
                  Ia_zi.mode          = posterior.mode(as.mcmc(post_zi[["Ia.obs"]])),
                  Ia_zi.mean          = mean(as.mcmc(post_zi[["Ia.obs"]])),
                  Ia_zi.low           = HPDinterval(as.mcmc(post_zi[["Ia.obs"]]))[1],
                  Ia_zi.high          = HPDinterval(as.mcmc(post_zi[["Ia.obs"]]))[2])

