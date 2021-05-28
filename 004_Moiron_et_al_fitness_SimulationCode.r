# Code for "The quantitative genetics of fitness in a wild bird population"
# Unpublished manuscript, doi: tba
# Moiron M, Charmantier A, Bouwhuis S

# The code provided here is sufficient to replicate the simulations presented in the above paper, AND
# replicates the code provided by de Villemereuil et al (2019) in the manuscript
# with the title "Little Adaptive Potential in a Threatened Passerine Bird" published in Current Biology.
# As such, please, use and cite the original source:
# de Villemereuil et al (2019) MS doi: 10.1016/j.cub.2019.01.072
# de Villemereuil et al (2019) code doi: 10.14324/000.ds.10065966


######################################################
# DATA SIMULATIONs
######################################################

# Loading packages
library(tidyverse)
library(MCMCglmm)
library(QGglmm)
library(pedantics)
library(beepr)

# Loading phenotypic data and pedigree
data <- read.table("fitness_data.txt", header=TRUE)
pedigree<- read.table("pedigree.txt",header=TRUE)

# Retrieving year of birth
years <- data[["year"]][match(pedigree[["animal"]], data[["chickID"]])]
table(years)

# Setting same year for all "founders"
years[is.na(years)] <- "Year0"

# Setting simulation parameters
K <- 100  # Number of simulations
N <- nrow(pedigree) # Number of individuals    
    
# Since the meaning of evolvability for binomial traits is less obvious, 
# we simulate a heritability of 0.1 for the Zero-inflated component of fitness, 
# which would correspond to an evolvability of 0.03.

mu_zi <- 1.7        #Intercept (ZI)
va_zi <-1.5378      # Additive genetic variance (ZI)
vyear_zi <-  0.138  # Year effect (ZI)  
vr_zi <- 1          # Residual variance (ZI)  

# Since heritability is not the best measure for additive genetic variance of lifetime fitnes, 
# we used an evolvability (or additive genetic variance of relative fitness) of 0.01 for the Poisson component

mu_pois <- 1.27    # Intercept (Poisson)
va_pois <- 0.01    # Additive genetic variance (Poisson)
vyear_pois <-0.12  # Year effect (Poisson)    
vr_pois <- 0.12    # Residual variance (Poisson)    

va_pois/(va_pois+vyear_pois +vr_pois) ##should correspond to a heritability of 0.04

# Simulating the data
phens <- matrix(0, nrow = N, ncol = K)

for (k in 1:K) {
    ## Poisson
    # Simulating the breeding values
    a_pois <- rbv(pedigree, va_pois)
    # Simulating the year effect
    y_eff_pois <- sapply(unique(years), function(string) {rnorm(1, 0, sqrt(vyear_pois))})
    # Latent trait for Poisson
    lat_pois <- mu_pois + 
                a_pois + 
                y_eff_pois[years] + 
                rnorm(N, 0, sqrt(vr_pois))
    lambda <- exp(lat_pois)
    
    ## Zero-inflation
    # Simulating the breeding values
    a_zi <- rbv(pedigree, va_zi)
    # Simulating the year effect
    y_eff_zi <- sapply(unique(years), function(string) {rnorm(1, 0, sqrt(vyear_zi))})
    # Latent trait for ZI
    lat_zi <- mu_zi + 
              a_zi + 
              y_eff_zi[years] + 
              rnorm(N, 0, sqrt(vr_zi))
    p_zi <- plogis(lat_zi)
    
    # Simulating the ZI-Poisson phenotype
    phens[ , k] <- ifelse(runif(N) < p_zi,
                          0,
                          rpois(N, lambda = lambda))
}

# Checking distribution of simulated data
hist(phens, breaks=35, xlab = "Phenotype",main = "")

# Analysing the data
colnames(pedigree) <- c("animal", "dam", "sire")

# Defining prior
prior <- list(R = list(V = diag(2), nu = 2,  fix = 2),
              G = list(G1 = list(V = diag(2)/2, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                       G2 = list(V =diag(2)/2, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000)))


# Creating the output dataframe
out <- data.frame(sim = 1:K,
                  va_pois.mode  = NA,
                  va_pois.mean= NA,
                  va_pois.low = NA,
                  va_pois.high= NA,
                  vp_pois.mode= NA,
                  vp_pois.mean= NA,
                  vp_pois.low = NA,
                  vp_pois.high= NA,
                  h2_pois.mode= NA,
                  h2_pois.mean= NA,
                  h2_pois.low = NA,
                  h2_pois.high= NA,
                  Ia_pois.mode= NA,
                  Ia_pois.mean= NA,
                  Ia_pois.low = NA,
                  Ia_pois.high= NA,
                  va_zi.mode  = NA,
                  va_zi.mean  = NA,
                  va_zi.low   = NA,
                  va_zi.high  = NA,
                  vp_zi.mode  = NA,
                  vp_zi.mean  = NA,
                  vp_zi.low   = NA,
                  vp_zi.high  = NA,
                  h2_zi.mode  = NA,
                  h2_zi.mean  = NA,
                  h2_zi.low   = NA,
                  h2_zi.high  = NA,
                  Ia_zi.mode  = NA,
                  Ia_zi.mean  = NA,
                  Ia_zi.low   = NA,
                  Ia_zi.high  = NA)

# Looping over the models
for(k in 1:K) {    
    # Creating a dataset
    data <- data.frame(phen = phens[ , k], 
                       animal = pedigree[["animal"]],
                       year = years)

    print(paste0("Performing model number ",k,"..."))
    
    # Running the model
    mod <- MCMCglmm(phen ~ trait - 1,
                   random = ~ idh(trait):animal +
                              idh(trait):year,
                   rcov = ~ idh(trait):units,
                   data = data, pedigree = pedigree, prior = prior, family = "zipoisson",
                   nitt = 101000, thin = 10, burnin = 1000, verbose = FALSE)

    # Computing h2 and Ia for ZI component
    df_zi <- data.frame(mu = as.vector(mod[["Sol"]][ , "traitzi_phen"]),
                        va = as.vector(mod[["VCV"]][ , "traitzi_phen.animal"]),
                        vp = rowSums(mod[["VCV"]][ , grepl("traitzi_phen", colnames(mod[["VCV"]]))]))
    post_zi <-  do.call("rbind",apply(df_zi, 1, function(row) {
        QGparams(mu = row["mu"], var.a = row["va"], var.p = row["vp"],
                 model = "binom1.logit", verbose = FALSE)
    }))
    post_zi[["Ia.obs"]] <- post_zi[["var.a.obs"]] / (post_zi[["mean.obs"]]^2)
  
    # Computing h2 and Ia for Poisson component
    df_pois <- data.frame(mu = as.vector(mod[["Sol"]][ , "traitphen"]),
                          va = as.vector(mod[["VCV"]][ , "traitphen.animal"]),
                          vp = rowSums(mod[["VCV"]][ , grepl("traitphen", colnames(mod[["VCV"]]))]))
    post_pois <-  do.call("rbind",apply(df_pois, 1, function(row) {
        QGparams(mu = row["mu"], var.a = row["va"], var.p = row["vp"],
                 model = "Poisson.log", verbose = FALSE)
    }))
    post_pois[["Ia.obs"]] <- post_pois[["var.a.obs"]] / (post_pois[["mean.obs"]]^2)
    
    # Saving model output
    out[k, ] <- data.frame(sim = k,
               va_pois.mode      = posterior.mode(as.mcmc(post_pois[["var.a.obs"]])),
               va_pois.mean      = mean(as.mcmc(post_pois[["var.a.obs"]])),
               va_pois.low       = HPDinterval(as.mcmc(post_pois[["var.a.obs"]]))[1],
               va_pois.high      = HPDinterval(as.mcmc(post_pois[["var.a.obs"]]))[2],
               vp_pois.mode      = posterior.mode(as.mcmc(post_pois[["var.obs"]])),
               vp_pois.mean      = mean(as.mcmc(post_pois[["var.obs"]])),
               vp_pois.low       = HPDinterval(as.mcmc(post_pois[["var.obs"]]))[1],
               vp_pois.high      = HPDinterval(as.mcmc(post_pois[["var.obs"]]))[2],
               h2_pois.mode      = posterior.mode(as.mcmc(post_pois[["h2.obs"]])),
               h2_pois.mean      = mean(as.mcmc(post_pois[["h2.obs"]])),
               h2_pois.low       = HPDinterval(as.mcmc(post_pois[["h2.obs"]]))[1],
               h2_pois.high      = HPDinterval(as.mcmc(post_pois[["h2.obs"]]))[2],
               Ia_pois.mode      = posterior.mode(as.mcmc(post_pois[["Ia.obs"]])),
               Ia_pois.mean      = mean(as.mcmc(post_pois[["Ia.obs"]])),
               Ia_pois.low       = HPDinterval(as.mcmc(post_pois[["Ia.obs"]]))[1],
               Ia_pois.high      = HPDinterval(as.mcmc(post_pois[["Ia.obs"]]))[2],
               va_zi.mode        = posterior.mode(as.mcmc(post_zi[["var.a.obs"]])),
               va_zi.mean        = mean(as.mcmc(post_zi[["var.a.obs"]])),
               va_zi.low         = HPDinterval(as.mcmc(post_zi[["var.a.obs"]]))[1],
               va_zi.high        = HPDinterval(as.mcmc(post_zi[["var.a.obs"]]))[2],
               vp_zi.mode        = posterior.mode(as.mcmc(post_zi[["var.obs"]])),
               vp_zi.mean        = mean(as.mcmc(post_zi[["var.obs"]])),
               vp_zi.low         = HPDinterval(as.mcmc(post_zi[["var.obs"]]))[1],
               vp_zi.high        = HPDinterval(as.mcmc(post_zi[["var.obs"]]))[2],
               h2_zi.mode        = posterior.mode(as.mcmc(post_zi[["h2.obs"]])),
               h2_zi.mean        = mean(as.mcmc(post_zi[["h2.obs"]])),
               h2_zi.low         = HPDinterval(as.mcmc(post_zi[["h2.obs"]]))[1],
               h2_zi.high        = HPDinterval(as.mcmc(post_zi[["h2.obs"]]))[2],
               Ia_zi.mode        = posterior.mode(as.mcmc(post_zi[["Ia.obs"]])),
               Ia_zi.mean        = mean(as.mcmc(post_zi[["Ia.obs"]])),
               Ia_zi.low         = HPDinterval(as.mcmc(post_zi[["Ia.obs"]]))[1],
               Ia_zi.high        = HPDinterval(as.mcmc(post_zi[["Ia.obs"]]))[2])
}
    
#save(out, file = "power_analysis_ziPoisson.Rdata")
#load("power_analysis_ziPoisson.Rdata")
    

######################################################
# Poisson component

    round(mean(out$Ia_pois.mean),3) #should be 0.01
    round(mean(out$Ia_pois.mode),3)
    
    round(mean(out$h2_pois.mean),3) #should be 0.04
    round(mean(out$h2_pois.mode),3)


######################################################
# ZI component
    
    round(mean(out$Ia_zi.mean),2) #should be 0.03
    round(mean(out$Ia_zi.mode),2)
    
    
    round(mean(out$h2_zi.mean),3) #should be 0.1
    round(mean(out$h2_zi.mode),3)
    
