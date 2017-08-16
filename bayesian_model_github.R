#################
# Bayesian inference on NESARC data using BAAHMS priors
# Citation: Goldstein ND, Burstyn I, Welles SL. Bayesian approaches to racial disparities in HIV risk estimation among men who have sex with men. Epidemiology. 2017 Mar;28(2):215-220.
# 12/10/14 -- Neal Goldstein
#################

### FUNCTIONS ###

library("rjags")
library("dclone") #parallel mcmc
library("boot") #bootstrapping

#boot strap prevalence of MSM from naive data
#called from boot with arguments dataframe, index of samples to pull, and model to run
bootPrev = function(data, index, indicator)
{
  newdata = data[index,]
  if (indicator=="identity")
  {
    return(sum(newdata$sex_identity_gay, na.rm=T)/length(na.omit(newdata$sex_identity_gay)))
  }
  else
  {
    return(sum(newdata$sex_men, na.rm=T)/length(na.omit(newdata$sex_men)))
  }
}
  

### READ DATA ###

load("NESARCw2.RData")


### COMPLETE CASE ###

#black
NESARCw2_partner = na.omit(NESARCw2_black[,c("sex_men","hiv_aids","relationship","std","abuse_sexual","recent_drug")])

#white
NESARCw2_partner = na.omit(NESARCw2_white[,c("sex_men","hiv_aids","relationship","std","abuse_sexual","recent_drug")])

#all
NESARCw2_partner = na.omit(NESARCw2_men[NESARCw2_men$black==1 | NESARCw2_men$white==1,c("sex_men","black","hiv_aids","relationship","std","abuse_sexual","recent_drug")])


### BUGS MODEL ###

#all men, latent confounding interaction model
bugs_model =
"model {
  for (i in 1:n) {
  
    #outcome model, log odds of hiv_aids given these predictors
    hiv_aids[i] ~ dbern(p_hiv_aids[i])
    logit(p_hiv_aids[i]) <- b0+b1*msm[i]+b2*black[i]+b3*msm[i]*black[i]+b4*relationship[i]+b5*std[i]+b6*abuse_sexual[i]+b7*recent_drug[i]+b8*U[i]
    
    #exposure models, log odds of true msm status given these predictors
    msm[i] ~ dbern(p_msm[i])
    logit(p_msm[i]) <- a0+a1*black[i]*a2*relationship[i]+a3*std[i]+a4*abuse_sexual[i]+a5*recent_drug[i]+a6*U[i]
    
    #measurement model, imputing the true msm status given the measurement error
    #allows for differential misclassification by HIV status as well as different estimates by race
    msm.star[i] ~ dbern(p_msm.star[i])
    p_msm.star[i] <- black[i]*(sn.black.msm.hivneg*msm[i]*(1-hiv_aids[i])+(1-msm[i])*(1-sp.black.msm.hivneg)*(1-hiv_aids[i]) + sn.black.msm.hivpos*msm[i]*(hiv_aids[i])+(1-msm[i])*(1-sp.black.msm.hivpos)*(hiv_aids[i])) + (1-black[i])*(sn.nonblack.msm.hivneg*msm[i]*(1-hiv_aids[i])+(1-msm[i])*(1-sp.nonblack.msm.hivneg)*(1-hiv_aids[i]) + sn.nonblack.msm.hivpos*msm[i]*(hiv_aids[i])+(1-msm[i])*(1-sp.nonblack.msm.hivpos)*(hiv_aids[i]))
    
    #prevalence models of potential confounders
    black[i] ~ dbern(p_black[i])
    logit(p_black[i]) <- prev.black

    relationship[i] ~ dbern(p_relationship[i])
    logit(p_relationship[i]) <- prev.relationship
    
    std[i] ~ dbern(p_std[i])
    logit(p_std[i]) <- prev.std
    
    abuse_sexual[i] ~ dbern(p_abuse_sexual[i])
    logit(p_abuse_sexual[i]) <- prev.abuse_sexual
    
    recent_drug[i] ~ dbern(p_recent_drug[i])
    logit(p_recent_drug[i]) <- prev.recent_drug
    
    #prevalence model of unknown counfounder
    U[i] ~ dbern(p_U[i])
    logit(p_U[i]) <- prev.U
  }
  
  #priors
  #for normal distribution, provide (mean, precision=(1/variance))
  #for beta distribution, provide (alpha, beta)
  
  #for independent priors use dnorm(0,1/(fixed variance, e.g. 10)
  #for hierarchical priors use dnorm(0,1/(random variance sample from inverse chi sq distribution); see: http://www.ncbi.nlm.nih.gov/pubmed/18226747
  #for prevalence priors, use independent prior dnorm(0,1/10)
  
  #instead of inverse chi2 (not in JAGS) use gamma (do not need to take recip), see: http://www.cs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf, and http://www.stat.ubc.ca/~gavin/WinBUGSdocs/WinBUGS%20lectures%20.pdf
  xi ~ dgamma(10,((log(6)/2)^2)) # hierarchical prior variance for outcome vars
  rho ~ dgamma(10,((log(6)/2)^2)) # hierarchical prior variance for exposure vars
  
  b0 ~ dnorm(0,xi)
  #b1 ~ dnorm(1.81,1/0.20) #informed prior from BAAHMS
  #b1 ~ dnorm(0,xi) #flat prior
  b1 ~ dnorm(1.81,xi) #informed prior from BAAHMS
  b2 ~ dnorm(0,xi)
  b3 ~ dnorm(0,xi)
  b4 ~ dnorm(0,xi)
  b5 ~ dnorm(0,xi)
  b6 ~ dnorm(0,xi)
  b7 ~ dnorm(0,xi)
  b8 ~ dnorm(0,xi) # log odds for relationship: U --> Y
  a0 ~ dnorm(0,rho)
  a1 ~ dnorm(0,rho)
  a2 ~ dnorm(0,rho)
  a3 ~ dnorm(0,rho)
  a4 ~ dnorm(0,rho)
  a5 ~ dnorm(0,rho)
  a6 ~ dnorm(0, rho) # log odds for relationship: U --> X
  prev.black ~ dnorm(0,1/10)
  prev.relationship ~ dnorm(0,1/10)
  prev.std ~ dnorm(0,1/10)
  prev.abuse_sexual ~ dnorm(0,1/10)
  prev.recent_drug ~ dnorm(0,1/10)
  prev.U ~ dnorm(0,1/10)

  #can have different SN/SP for black,nonblack men although keep it same right now
  sn.black.msm.hivneg ~ dbeta(45,6) #partner, add beta(1,1)
  sp.black.msm.hivneg ~ dbeta(518,20) #partner, add beta(1,1)
  sn.black.msm.hivpos ~ dbeta(15,2) #partner, add beta(1,1)
  sp.black.msm.hivpos ~ dbeta(17,1) #partner, add beta(1,1)
  sn.nonblack.msm.hivneg ~ dbeta(45,6) #partner, add beta(1,1)
  sp.nonblack.msm.hivneg ~ dbeta(518,20) #partner, add beta(1,1)
  sn.nonblack.msm.hivpos ~ dbeta(15,2) #partner, add beta(1,1)
  sp.nonblack.msm.hivpos ~ dbeta(17,1) #partner, add beta(1,1)
}"


#all men, latent confounding model
bugs_model =
  "model {
    for (i in 1:n) {
    
      #outcome model, log odds of hiv_aids given these predictors
      hiv_aids[i] ~ dbern(p_hiv_aids[i])
      logit(p_hiv_aids[i]) <- b0+b1*msm[i]+b2*black[i]+b4*relationship[i]+b5*std[i]+b6*abuse_sexual[i]+b7*recent_drug[i]+b8*U[i]
      
      #exposure models, log odds of true msm status given these predictors
      msm[i] ~ dbern(p_msm[i])
      logit(p_msm[i]) <- a0+a1*black[i]*a2*relationship[i]+a3*std[i]+a4*abuse_sexual[i]+a5*recent_drug[i]+a6*U[i]
      
      #measurement model, imputing the true msm status given the measurement error
      #allows for differential misclassification by HIV status as well as different estimates by race
      msm.star[i] ~ dbern(p_msm.star[i])
      p_msm.star[i] <- black[i]*(sn.black.msm.hivneg*msm[i]*(1-hiv_aids[i])+(1-msm[i])*(1-sp.black.msm.hivneg)*(1-hiv_aids[i]) + sn.black.msm.hivpos*msm[i]*(hiv_aids[i])+(1-msm[i])*(1-sp.black.msm.hivpos)*(hiv_aids[i])) + (1-black[i])*(sn.nonblack.msm.hivneg*msm[i]*(1-hiv_aids[i])+(1-msm[i])*(1-sp.nonblack.msm.hivneg)*(1-hiv_aids[i]) + sn.nonblack.msm.hivpos*msm[i]*(hiv_aids[i])+(1-msm[i])*(1-sp.nonblack.msm.hivpos)*(hiv_aids[i]))
      
      #prevalence models of potential confounders
      black[i] ~ dbern(p_black[i])
      logit(p_black[i]) <- prev.black
      
      relationship[i] ~ dbern(p_relationship[i])
      logit(p_relationship[i]) <- prev.relationship
      
      std[i] ~ dbern(p_std[i])
      logit(p_std[i]) <- prev.std
      
      abuse_sexual[i] ~ dbern(p_abuse_sexual[i])
      logit(p_abuse_sexual[i]) <- prev.abuse_sexual
      
      recent_drug[i] ~ dbern(p_recent_drug[i])
      logit(p_recent_drug[i]) <- prev.recent_drug
      
      #prevalence model of unknown counfounder
      U[i] ~ dbern(p_U[i])
      logit(p_U[i]) <- prev.U
    }

    #priors
    #for normal distribution, provide (mean, precision=(1/variance))
    #for beta distribution, provide (alpha, beta)
    
    #for independent priors use dnorm(0,1/(fixed variance, e.g. 10)
    #for hierarchical priors use dnorm(0,1/(random variance sample from inverse chi sq distribution); see: http://www.ncbi.nlm.nih.gov/pubmed/18226747
    #for prevalence priors, use independent prior dnorm(0,1/10)
    
    #instead of inverse chi2 (not in JAGS) use gamma (do not need to take recip), see: http://www.cs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf, and http://www.stat.ubc.ca/~gavin/WinBUGSdocs/WinBUGS%20lectures%20.pdf
    xi ~ dgamma(10,((log(6)/2)^2)) # hierarchical prior variance for outcome vars
    rho ~ dgamma(10,((log(6)/2)^2)) # hierarchical prior variance for exposure vars
    
    b0 ~ dnorm(0,xi)
    #b1 ~ dnorm(1.81,1/0.20) #informed prior from BAAHMS
    #b1 ~ dnorm(0,xi) #flat prior
    b1 ~ dnorm(1.81,xi) #informed prior from BAAHMS
    b2 ~ dnorm(0,xi)
    b4 ~ dnorm(0,xi)
    b5 ~ dnorm(0,xi)
    b6 ~ dnorm(0,xi)
    b7 ~ dnorm(0,xi)
    b8 ~ dnorm(0,xi) # log odds for relationship: U --> Y
    a0 ~ dnorm(0,rho)
    a1 ~ dnorm(0,rho)
    a2 ~ dnorm(0,rho)
    a3 ~ dnorm(0,rho)
    a4 ~ dnorm(0,rho)
    a5 ~ dnorm(0,rho)
    a6 ~ dnorm(0, rho) # log odds for relationship: U --> X
    prev.black ~ dnorm(0,1/10)
    prev.relationship ~ dnorm(0,1/10)
    prev.std ~ dnorm(0,1/10)
    prev.abuse_sexual ~ dnorm(0,1/10)
    prev.recent_drug ~ dnorm(0,1/10)
    prev.U ~ dnorm(0,1/10)
    
    #can have different SN/SP for black,nonblack men although keep it same right now
    sn.black.msm.hivneg ~ dbeta(45,6) #partner, add beta(1,1)
    sp.black.msm.hivneg ~ dbeta(518,20) #partner, add beta(1,1)
    sn.black.msm.hivpos ~ dbeta(15,2) #partner, add beta(1,1)
    sp.black.msm.hivpos ~ dbeta(17,1) #partner, add beta(1,1)
    sn.nonblack.msm.hivneg ~ dbeta(45,6) #partner, add beta(1,1)
    sp.nonblack.msm.hivneg ~ dbeta(518,20) #partner, add beta(1,1)
    sn.nonblack.msm.hivpos ~ dbeta(15,2) #partner, add beta(1,1)
    sp.nonblack.msm.hivpos ~ dbeta(17,1) #partner, add beta(1,1)
}"

#black MSM, latent confounding stratified model
bugs_model =
"model {
  for (i in 1:n) {
    
    #outcome model, log odds of hiv_aids given these predictors
    hiv_aids[i] ~ dbern(p_hiv_aids[i])
    logit(p_hiv_aids[i]) <- b0+b1*msm[i]+b2*relationship[i]+b3*std[i]+b4*abuse_sexual[i]+b5*recent_drug[i]+b6*U[i]
    
    #exposure models, log odds of true msm status given these predictors
    msm[i] ~ dbern(p_msm[i])
    logit(p_msm[i]) <- a0+a1*relationship[i]+a2*std[i]+a3*abuse_sexual[i]+a4*recent_drug[i]+a5*U[i]
    
    #measurement model, imputing the true msm status given the measurement error
    msm.star[i] ~ dbern(p_msm.star[i])
    p_msm.star[i] <- sn.msm.hivneg*msm[i]*(1-hiv_aids[i])+(1-msm[i])*(1-sp.msm.hivneg)*(1-hiv_aids[i]) + sn.msm.hivpos*msm[i]*(hiv_aids[i])+(1-msm[i])*(1-sp.msm.hivpos)*(hiv_aids[i])
    
    #prevalence models of potential confounders
    relationship[i] ~ dbern(p_relationship[i])
    logit(p_relationship[i]) <- prev.relationship
    
    std[i] ~ dbern(p_std[i])
    logit(p_std[i]) <- prev.std
    
    abuse_sexual[i] ~ dbern(p_abuse_sexual[i])
    logit(p_abuse_sexual[i]) <- prev.abuse_sexual
    
    recent_drug[i] ~ dbern(p_recent_drug[i])
    logit(p_recent_drug[i]) <- prev.recent_drug
    
    #prevalence model of unknown counfounder
    U[i] ~ dbern(p_U[i])
    logit(p_U[i]) <- prev.U
  }
  
  #priors
  #for normal distribution, provide (mean, precision=(1/variance))
  #for beta distribution, provide (alpha, beta)
  
  #for independent priors use dnorm(0,1/(fixed variance, e.g. 10)
  #for hierarchical priors use dnorm(0,1/(random variance sample from inverse chi sq distribution); see: http://www.ncbi.nlm.nih.gov/pubmed/18226747
  #for prevalence priors, use independent prior dnorm(0,1/10)

  #instead of inverse chi2 (not in JAGS) use gamma (do not need to take recip), see: http://www.cs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf, and http://www.stat.ubc.ca/~gavin/WinBUGSdocs/WinBUGS%20lectures%20.pdf
  xi ~ dgamma(10,((log(6)/2)^2)) # hierarchical prior variance for outcome vars
  rho ~ dgamma(10,((log(6)/2)^2)) # hierarchical prior variance for exposure vars

  b0 ~ dnorm(0,xi)
  #b1 ~ dnorm(1.81,1/0.20) #informed prior from BAAHMS
  #b1 ~ dnorm(0,xi) #flat prior
  b1 ~ dnorm(1.81,xi) #informed prior from BAAHMS
  b2 ~ dnorm(0,xi)
  b3 ~ dnorm(0,xi)
  b4 ~ dnorm(0,xi)
  b5 ~ dnorm(0,xi)
  b6 ~ dnorm(0,xi) # log odds for relationship: U --> Y
  a0 ~ dnorm(0,rho)
  a1 ~ dnorm(0,rho)
  a2 ~ dnorm(0,rho)
  a3 ~ dnorm(0,rho)
  a4 ~ dnorm(0,rho)
  a5 ~ dnorm(0,rho) # log odds for relationship: U --> X
  prev.relationship ~ dnorm(0,1/10)
  prev.std ~ dnorm(0,1/10)
  prev.abuse_sexual ~ dnorm(0,1/10)
  prev.recent_drug ~ dnorm(0,1/10)
  prev.U ~ dnorm(0,1/10)
  sn.msm.hivneg ~ dbeta(45,6) #partner, add beta(1,1)
  sp.msm.hivneg ~ dbeta(518,20) #partner, add beta(1,1)
  sn.msm.hivpos ~ dbeta(15,2) #partner, add beta(1,1)
  sp.msm.hivpos ~ dbeta(17,1) #partner, add beta(1,1)
}"

#white MSM, latent confounding stratified model
bugs_model =
"model {
  for (i in 1:n) {
  
    #outcome model, log odds of hiv_aids given these predictors
    hiv_aids[i] ~ dbern(p_hiv_aids[i])
    logit(p_hiv_aids[i]) <- b0+b1*msm[i]+b2*relationship[i]+b3*std[i]+b4*abuse_sexual[i]+b5*recent_drug[i]+b6*U[i]
    
    #exposure models, log odds of true msm status given these predictors
    msm[i] ~ dbern(p_msm[i])
    logit(p_msm[i]) <- a0+a1*relationship[i]+a2*std[i]+a3*abuse_sexual[i]+a4*recent_drug[i]+a5*U[i]
    
    #measurement model, imputing the true msm status given the measurement error
    msm.star[i] ~ dbern(p_msm.star[i])
    p_msm.star[i] <- sn.msm.hivneg*msm[i]*(1-hiv_aids[i])+(1-msm[i])*(1-sp.msm.hivneg)*(1-hiv_aids[i]) + sn.msm.hivpos*msm[i]*(hiv_aids[i])+(1-msm[i])*(1-sp.msm.hivpos)*(hiv_aids[i])
    
    #prevalence models of potential confounders
    relationship[i] ~ dbern(p_relationship[i])
    logit(p_relationship[i]) <- prev.relationship
    
    std[i] ~ dbern(p_std[i])
    logit(p_std[i]) <- prev.std
    
    abuse_sexual[i] ~ dbern(p_abuse_sexual[i])
    logit(p_abuse_sexual[i]) <- prev.abuse_sexual
    
    recent_drug[i] ~ dbern(p_recent_drug[i])
    logit(p_recent_drug[i]) <- prev.recent_drug
  
    #prevalence model of unknown counfounder
    U[i] ~ dbern(p_U[i])
    logit(p_U[i]) <- prev.U
  }
  
  #priors
  #for normal distribution, provide (mean, precision=(1/variance))
  #for beta distribution, provide (alpha, beta)
  
  #for independent priors use dnorm(0,1/(fixed variance, e.g. 10)
  #for hierarchical priors use dnorm(0,1/(random variance sample from inverse chi sq distribution); see: http://www.ncbi.nlm.nih.gov/pubmed/18226747
  #for prevalence priors, use independent prior dnorm(0,1/10)

  #instead of inverse chi2 (not in JAGS) use gamma (do not need to take recip), see: http://www.cs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf, and http://www.stat.ubc.ca/~gavin/WinBUGSdocs/WinBUGS%20lectures%20.pdf
  xi ~ dgamma(10,((log(6)/2)^2)) # hierarchical prior variance for outcome vars
  rho ~ dgamma(10,((log(6)/2)^2)) # hierarchical prior variance for exposure vars

  b0 ~ dnorm(0,xi)
  #b1 ~ dnorm(1.81,1/0.20) #informed prior from BAAHMS
  #b1 ~ dnorm(0,xi) #flat prior
  b1 ~ dnorm(1.81,xi) #informed prior from BAAHMS
  b2 ~ dnorm(0,xi)
  b3 ~ dnorm(0,xi)
  b4 ~ dnorm(0,xi)
  b5 ~ dnorm(0,xi)
  b6 ~ dnorm(0,xi) # log odds for relationship: U --> Y
  a0 ~ dnorm(0,rho)
  a1 ~ dnorm(0,rho)
  a2 ~ dnorm(0,rho)
  a3 ~ dnorm(0,rho)
  a4 ~ dnorm(0,rho)
  a5 ~ dnorm(0,rho) # log odds for relationship: U --> X
  prev.relationship ~ dnorm(0,1/10)
  prev.std ~ dnorm(0,1/10)
  prev.abuse_sexual ~ dnorm(0,1/10)
  prev.recent_drug ~ dnorm(0,1/10)
  prev.U ~ dnorm(0,1/10)
  sn.msm.hivneg ~ dbeta(45,6) #partner, add beta(1,1)
  #sn.msm.hivneg ~ dbeta(48,3) #increase to ~95% (actually 95.9%, move 3 from no to yes)
  sp.msm.hivneg ~ dbeta(518,20) #partner, add beta(1,1)
  sn.msm.hivpos ~ dbeta(15,2) #partner, add beta(1,1)
  #sn.msm.hivpos ~ dbeta(15,2) #increase to 95% (actually 93.3%, can't increase)
  sp.msm.hivpos ~ dbeta(17,1) #partner, add beta(1,1)
}"

#black MSM, misclassification alone stratified model
bugs_model =
"model {
  for (i in 1:n) {
  
    #outcome model, log odds of hiv_aids given these predictors
    hiv_aids[i] ~ dbern(p_hiv_aids[i])
    logit(p_hiv_aids[i]) <- b0+b1*msm[i]+b2*relationship[i]+b3*std[i]+b4*abuse_sexual[i]+b5*recent_drug[i]
    
    #exposure models, log odds of true msm status given these predictors
    msm[i] ~ dbern(p_msm[i])
    logit(p_msm[i]) <- a0+a1*relationship[i]+a2*std[i]+a3*abuse_sexual[i]+a4*recent_drug[i]
    
    #measurement model, imputing the true msm status given the measurement error
    msm.star[i] ~ dbern(p_msm.star[i])
    p_msm.star[i] <- sn.msm.hivneg*msm[i]*(1-hiv_aids[i])+(1-msm[i])*(1-sp.msm.hivneg)*(1-hiv_aids[i]) + sn.msm.hivpos*msm[i]*(hiv_aids[i])+(1-msm[i])*(1-sp.msm.hivpos)*(hiv_aids[i])
    
    #prevalence models of potential confounders
    relationship[i] ~ dbern(p_relationship[i])
    logit(p_relationship[i]) <- prev.relationship
    
    std[i] ~ dbern(p_std[i])
    logit(p_std[i]) <- prev.std
    
    abuse_sexual[i] ~ dbern(p_abuse_sexual[i])
    logit(p_abuse_sexual[i]) <- prev.abuse_sexual
    
    recent_drug[i] ~ dbern(p_recent_drug[i])
    logit(p_recent_drug[i]) <- prev.recent_drug
  }
  
  #priors
  #for normal distribution, provide (mean, precision=(1/variance))
  #for beta distribution, provide (alpha, beta)

  b0 ~ dnorm(0,1/10)
  b1 ~ dnorm(1.81,1/0.20) #informed prior from BAAHMS
  b1 ~ dnorm(0,1/10)
  b2 ~ dnorm(0,1/10)
  b3 ~ dnorm(0,1/10)
  b4 ~ dnorm(0,1/10)
  b5 ~ dnorm(0,1/10)
  a0 ~ dnorm(0,1/10)
  a1 ~ dnorm(0,1/10)
  a2 ~ dnorm(0,1/10)
  a3 ~ dnorm(0,1/10)
  a4 ~ dnorm(0,1/10)
  prev.relationship ~ dnorm(0,1/10)
  prev.std ~ dnorm(0,1/10)
  prev.abuse_sexual ~ dnorm(0,1/10)
  prev.recent_drug ~ dnorm(0,1/10)
  sn.msm.hivneg ~ dbeta(45,6) #partner, add beta(1,1)
  sp.msm.hivneg ~ dbeta(518,20) #partner, add beta(1,1)
  sn.msm.hivpos ~ dbeta(15,2) #partner, add beta(1,1)
  sp.msm.hivpos ~ dbeta(17,1) #partner, add beta(1,1)
}"

#white MSM, misclassification alone stratified model
bugs_model =
"model {
  for (i in 1:n) {
    
    #outcome model, log odds of hiv_aids given these predictors
    hiv_aids[i] ~ dbern(p_hiv_aids[i])
    logit(p_hiv_aids[i]) <- b0+b1*msm[i]+b2*relationship[i]+b3*std[i]+b4*abuse_sexual[i]+b5*recent_drug[i]

    #exposure models, log odds of true msm status given these predictors
    msm[i] ~ dbern(p_msm[i])
    logit(p_msm[i]) <- a0+a1*relationship[i]+a2*std[i]+a3*abuse_sexual[i]+a4*recent_drug[i]

    #measurement model, imputing the true msm status given the measurement error
    msm.star[i] ~ dbern(p_msm.star[i])
    p_msm.star[i] <- sn.msm.hivneg*msm[i]*(1-hiv_aids[i])+(1-msm[i])*(1-sp.msm.hivneg)*(1-hiv_aids[i]) + sn.msm.hivpos*msm[i]*(hiv_aids[i])+(1-msm[i])*(1-sp.msm.hivpos)*(hiv_aids[i])

    #prevalence models of potential confounders
    relationship[i] ~ dbern(p_relationship[i])
    logit(p_relationship[i]) <- prev.relationship

    std[i] ~ dbern(p_std[i])
    logit(p_std[i]) <- prev.std

    abuse_sexual[i] ~ dbern(p_abuse_sexual[i])
    logit(p_abuse_sexual[i]) <- prev.abuse_sexual

    recent_drug[i] ~ dbern(p_recent_drug[i])
    logit(p_recent_drug[i]) <- prev.recent_drug

  }

  #priors
  #for normal distribution, provide (mean, precision=(1/variance))
  #for beta distribution, provide (alpha, beta)

  b0 ~ dnorm(0,1/10)
  b1 ~ dnorm(1.81,1/0.20) #variance inflated by 0%
  #b1 ~ dnorm(1.81,1/0.25) #variance inflated by 25%
  #b1 ~ dnorm(1.81,1/0.30) #variance inflated by 50%
  #b1 ~ dnorm(1.81,1/0.35) #variance inflated by 75%
  #b1 ~ dnorm(1.81,1/0.40) #variance inflated by 100%
  b2 ~ dnorm(0,1/10)
  b3 ~ dnorm(0,1/10)
  b4 ~ dnorm(0,1/10)
  b5 ~ dnorm(0,1/10)
  a0 ~ dnorm(0,1/10)
  a1 ~ dnorm(0,1/10)
  a2 ~ dnorm(0,1/10)
  a3 ~ dnorm(0,1/10)
  a4 ~ dnorm(0,1/10)
  prev.relationship ~ dnorm(0,1/10)
  prev.std ~ dnorm(0,1/10)
  prev.abuse_sexual ~ dnorm(0,1/10)
  prev.recent_drug ~ dnorm(0,1/10)
  sn.msm.hivneg ~ dbeta(45,6) #partner, add beta(1,1)
  #sn.msm.hivneg ~ dbeta(48,3) #increase to ~95% (actually 95.9%, move 3 from no to yes)
  sp.msm.hivneg ~ dbeta(518,20) #partner, add beta(1,1)
  sn.msm.hivpos ~ dbeta(15,2) #partner, add beta(1,1)
  #sn.msm.hivpos ~ dbeta(15,2) #increase to 95% (actually 93.3%, can't increase)
  sp.msm.hivpos ~ dbeta(17,1) #partner, add beta(1,1)
}"


### BAYESIAN SAMPLING ###

#write bugs model to temp file for JAGS
writeLines(bugs_model, file("bugs_model.txt"))

#initialize cluster for two parallel chains
jags_cluster = makeCluster(2)

#initialize model
parJagsModel(jags_cluster, name="res", file="bugs_model.txt",
                   data = list('msm.star' = NESARCw2_partner$sex_men,
                               'hiv_aids' = NESARCw2_partner$hiv_aids,
                               'black' = NESARCw2_partner$black,                   
                               'relationship' = NESARCw2_partner$relationship,
                               'std' = NESARCw2_partner$std,
                               'abuse_sexual' = NESARCw2_partner$abuse_sexual,
                               'recent_drug' = NESARCw2_partner$recent_drug,
                               'n' = nrow(NESARCw2_partner)),
                   n.chains = 2,
                   n.adapt = 100)

#sample from the posterior distribution for stratified model
#jags_samples = parCodaSamples(jags_cluster, "res", variable.names=c("b1","b6","sn.msm.hivneg","sp.msm.hivneg","sn.msm.hivpos","sp.msm.hivpos","msm","U"), n.iter=10000)

#sample from the posterior distribution for all men model
jags_samples = parCodaSamples(jags_cluster, "res", variable.names=c("b1","b2","b8","sn.msm.hivneg","sp.msm.hivneg","sn.msm.hivpos","sp.msm.hivpos","msm","U"), n.iter=10000)

#sample from the posterior distribution for interaction model
#jags_samples = parCodaSamples(jags_cluster, "res", variable.names=c("b1","b2","b3","b8","sn.msm.hivneg","sp.msm.hivneg","sn.msm.hivpos","sp.msm.hivpos","msm","U"), n.iter=10000)

#clean up
stopCluster(jags_cluster)
rm(jags_cluster)
file.remove("bugs_model.txt")


### SAVE POSTERIORS ###

#black partner, MSM corrected
#save.image("Posterior distributions/Black partner.RData")
#load("Posterior distributions/Black partner.RData")

#black partner, unknown confounder, flat beta1
#save.image("Posterior distributions/Black partner latent flat.RData")
#load("Posterior distributions/Black partner latent flat.RData")

#black partner, unknown confounder, informed beta1
save.image("Posterior distributions/Black partner latent informed.RData")
load("Posterior distributions/Black partner latent informed.RData")

#white partner, MSM corrected
#save.image("Posterior distributions/Nonblack partner.RData")
#load("Posterior distributions/Nonblack partner.RData")

#white partner, unknown confounder, informed beta1
save.image("Posterior distributions/Nonblack partner latent informed.RData")
load("Posterior distributions/Nonblack partner latent informed.RData")

#all men partner, unknown confounder, informed beta1
save.image("Posterior distributions/All partner latent informed.RData")
load("Posterior distributions/All partner latent informed.RData")

#all men (interaction) partner, unknown confounder, informed beta1
save.image("Posterior distributions/All partner latent informed interaction.RData")
load("Posterior distributions/All partner latent informed interaction.RData")


### BAYESIAN INFERENCE ###

#check for convergence
plot(jags_samples[,"b1"]) #good mixing and shape of distribution
gelman.plot(jags_samples[,"b1"]) #no discernable difference
gelman.diag(jags_samples) #<1.1

#statistics for stratified model fit, discard first 1000 observations for burn in
#b1 = msm
#b6 = U
summary(window(jags_samples[,c("b1","b6")], start=1000))

#statistics for all men model fit, discard first 1000 observations for burn in
#b1 = msm
#b2 = black
#b3 = msm*black
#b8 = U
summary(window(jags_samples[,c("b1","b2","b8")], start=1000))

#statistics for interaction model fit, discard first 1000 observations for burn in
#b1 = msm
#b2 = black
#b3 = msm*black
#b8 = U
summary(window(jags_samples[,c("b1","b2","b3","b8")], start=1000))

#wald test for heterogeneity of effects; see PMID: 24022386

#estimates and standard errors
black = 1.737
black_se = 0.5255
white = 2.438
white_se = 0.5080

#compute test statistic (Z-score)
zscore = (white-black)/(sqrt((black_se^2) + (white_se^2)))

#check against normal distribution
2*pnorm(-abs(zscore))


### PREVALENCE of U ###

#stratified model

#obtain MSM status for each individual, for each simulation, in each chain
#each row represents a dataset, with each column being an individual
prevU_chain1 = as.data.frame(window(jags_samples[[1]], start=1000))
prevU_chain2 = as.data.frame(window(jags_samples[[2]], start=1000))

#drop unnecessary vars
prevU_chain1$b1 = NULL
prevU_chain1$b6 = NULL
prevU_chain1$sn.msm.hivneg = NULL
prevU_chain1$sn.msm.hivpos = NULL
prevU_chain1$sp.msm.hivneg = NULL
prevU_chain1$sp.msm.hivpos = NULL
prevU_chain2$b1 = NULL
prevU_chain2$b6 = NULL
prevU_chain2$sn.msm.hivneg = NULL
prevU_chain2$sn.msm.hivpos = NULL
prevU_chain2$sp.msm.hivneg = NULL
prevU_chain2$sp.msm.hivpos = NULL

#merge chains
prevU = rbind(prevU_chain1, prevU_chain2)

#separate U and MSM vars
tmpU = prevU[,1:(ncol(prevU)/2)]
tmpMSM = prevU[,((ncol(prevU)/2)+1):(ncol(prevU))]

#prev U | MSM
prevUMSM = rowSums(tmpU==1 & tmpMSM==1)/rowSums(tmpMSM)

#prev U | not MSM  
prevUnotMSM = rowSums(tmpU==1 & tmpMSM==0)/rowSums(!tmpMSM)

rm(tmpU,tmpMSM)

#summary of prev U | MSM
mean(prevUMSM)
quantile(prevUMSM,probs=c(0.025,0.975))

#summary of prev U | not MSM
mean(prevUnotMSM)
quantile(prevUnotMSM,probs=c(0.025,0.975))


