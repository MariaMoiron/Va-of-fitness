# Code for "The quantitative genetics of fitness in a wild bird population"
# Unpublished manuscript, doi: tba
# Moiron M, Charmantier A, Bouwhuis S

# The code provided here is sufficient to replicate the simulations presented in the above paper


######################################################
# DATA ANALYSIS OF ANNUAL REPRODUCTIVE SUCCESS DATA
######################################################

# Loading packages
library(tidyverse)
library(MCMCglmm)
library(QGglmm)
library(pedantics)
library(beepr)

# Loading phenotypic data
data <- read.table("fitness_data.txt", header=TRUE)

# Formating response varibale
data$ARS<-as.numeric(Data$ARS)  #annual reproductive success data

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
priorARS<- list(R = list(V = 1*0.03, nu = 1),
                      G = list(G1 = list(V =1*0.03, nu = 1, alpha.mu =0, alpha.V = 1000),
                               G2 = list(V = 1*0.1, nu = 1, alpha.mu =0, alpha.V = 1000),
                               G3 = list(V = 1*0.5, nu = 1, alpha.mu =0, alpha.V = 1000)))

# Running model
mod<- MCMCglmm(fitness ~ 1 +age,
               random = ~ animal +ID.PE +YEAR,
               rcov = ~ units,
               ginverse = list(animal = my_inverse),
               data = Data,
               prior =priorARS,
               family = "poisson",
               nitt = NITT, thin = THIN, burnin = BURN)

beep(1)
#save(mod, file = "analysis_ARS.rda")
#load("analysis_ARS.rda")

summary(mod)

# Assesing convergence 
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)
autocorr.diag(mod$VCV)

# Random effects in latent scale
round(posterior.mode(mod$VCV),3)
round(HPDinterval(mod$VCV), 3)
plot(mod$VCV)

# Computing h2 and Ia
df_pois <- data.frame(mu = as.vector(mod[["Sol"]][ , "(Intercept)"]),
                      va = as.vector(mod[["VCV"]][ , "animal"]),
                      vp = rowSums(mod[["VCV"]]))

post_pois <-  do.call("rbind",apply(df_pois, 1, function(row) {
  QGparams(mu = row["mu"], var.a = row["va"], var.p = row["vp"],
           model = "Poisson.log", verbose = FALSE)
}))

post_pois[["Ia.obs"]] <- post_pois[["var.a.obs"]] / (post_pois[["mean.obs"]]^2)

# Saving output
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
                  Ia_pois.high        = HPDinterval(as.mcmc(post_pois[["Ia.obs"]]))[2]
)