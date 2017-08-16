#################
# NESARC wave 2 analysis
# Citation: Goldstein ND, Burstyn I, Welles SL. Bayesian approaches to racial disparities in HIV risk estimation among men who have sex with men. Epidemiology. 2017 Mar;28(2):215-220.
# 7/7/14 -- Neal Goldstein
#################


### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(car) #vif
library(survey) #complex survey sampling
library(obsSens) #unknown confounder sensitivity analysis

### READ DATA ###

load("NESARCw2.RData")


### COMPLEX SURVEY ESTIMATES  ###

#some strata may have one sampled unit: http://r-survey.r-forge.r-project.org/survey/example-twostage.html
options(survey.lonely.psu="remove")

#specify survey design
NESARCw2_complex = svydesign(id=~psu, strata=~stratum, weights=~weight, data=NESARCw2_men)


### ANALYSIS: DESCRIPTIVES ###

#counts black
CrossTable(NESARCw2_black$sex_men)
CrossTable(NESARCw2_black$hiv_aids)
CrossTable(NESARCw2_black$relationship)
CrossTable(NESARCw2_black$std)
CrossTable(NESARCw2_black$abuse_sexual)
CrossTable(NESARCw2_black$recent_drug)

CrossTable(NESARCw2_black$hiv_aids, NESARCw2_black$sex_men, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_black$relationship, NESARCw2_black$sex_men, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_black$std, NESARCw2_black$sex_men, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_black$abuse_sexual, NESARCw2_black$sex_men, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_black$recent_drug, NESARCw2_black$sex_men, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

#proportions, black
svymean(~factor(black), NESARCw2_complex, na.rm=T)
svymean(~factor(sex_men), subset(NESARCw2_complex,black==1), na.rm=T)
svymean(~factor(hiv_aids), subset(NESARCw2_complex,black==1), na.rm=T)
svymean(~factor(hiv_aids), subset(NESARCw2_complex,black==1 & sex_men==1), na.rm=T)
svymean(~factor(hiv_aids), subset(NESARCw2_complex,black==1 & sex_men==0), na.rm=T)
svymean(~factor(relationship), subset(NESARCw2_complex,black==1), na.rm=T)
svymean(~factor(relationship), subset(NESARCw2_complex,black==1 & sex_men==1), na.rm=T)
svymean(~factor(relationship), subset(NESARCw2_complex,black==1 & sex_men==0), na.rm=T)
svymean(~factor(std), subset(NESARCw2_complex,black==1), na.rm=T)
svymean(~factor(std), subset(NESARCw2_complex,black==1 & sex_men==1), na.rm=T)
svymean(~factor(std), subset(NESARCw2_complex,black==1 & sex_men==0), na.rm=T)
svymean(~factor(abuse_sexual), subset(NESARCw2_complex,black==1), na.rm=T)
svymean(~factor(abuse_sexual), subset(NESARCw2_complex,black==1 & sex_men==1), na.rm=T)
svymean(~factor(abuse_sexual), subset(NESARCw2_complex,black==1 & sex_men==0), na.rm=T)
svymean(~factor(recent_drug), subset(NESARCw2_complex,black==1), na.rm=T)
svymean(~factor(recent_drug), subset(NESARCw2_complex,black==1 & sex_men==1), na.rm=T)
svymean(~factor(recent_drug), subset(NESARCw2_complex,black==1 & sex_men==0), na.rm=T)

#counts white
CrossTable(NESARCw2_white$sex_men)
CrossTable(NESARCw2_white$hiv_aids)
CrossTable(NESARCw2_white$relationship)
CrossTable(NESARCw2_white$std)
CrossTable(NESARCw2_white$abuse_sexual)
CrossTable(NESARCw2_white$recent_drug)

CrossTable(NESARCw2_white$hiv_aids, NESARCw2_white$sex_men, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_white$relationship, NESARCw2_white$sex_men, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_white$std, NESARCw2_white$sex_men, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_white$abuse_sexual, NESARCw2_white$sex_men, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_white$recent_drug, NESARCw2_white$sex_men, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

#proportions, black
svymean(~factor(white), NESARCw2_complex, na.rm=T)
svymean(~factor(sex_men), subset(NESARCw2_complex,white==1), na.rm=T)
svymean(~factor(hiv_aids), subset(NESARCw2_complex,white==1), na.rm=T)
svymean(~factor(hiv_aids), subset(NESARCw2_complex,white==1 & sex_men==1), na.rm=T)
svymean(~factor(hiv_aids), subset(NESARCw2_complex,white==1 & sex_men==0), na.rm=T)
svymean(~factor(relationship), subset(NESARCw2_complex,white==1), na.rm=T)
svymean(~factor(relationship), subset(NESARCw2_complex,white==1 & sex_men==1), na.rm=T)
svymean(~factor(relationship), subset(NESARCw2_complex,white==1 & sex_men==0), na.rm=T)
svymean(~factor(std), subset(NESARCw2_complex,white==1), na.rm=T)
svymean(~factor(std), subset(NESARCw2_complex,white==1 & sex_men==1), na.rm=T)
svymean(~factor(std), subset(NESARCw2_complex,white==1 & sex_men==0), na.rm=T)
svymean(~factor(abuse_sexual), subset(NESARCw2_complex,white==1), na.rm=T)
svymean(~factor(abuse_sexual), subset(NESARCw2_complex,white==1 & sex_men==1), na.rm=T)
svymean(~factor(abuse_sexual), subset(NESARCw2_complex,white==1 & sex_men==0), na.rm=T)
svymean(~factor(recent_drug), subset(NESARCw2_complex,white==1), na.rm=T)
svymean(~factor(recent_drug), subset(NESARCw2_complex,white==1 & sex_men==1), na.rm=T)
svymean(~factor(recent_drug), subset(NESARCw2_complex,white==1 & sex_men==0), na.rm=T)

#msm
CrossTable(NESARCw2_men$sex_men[NESARCw2_men$black==1 | NESARCw2_men$white==1])
svymean(~factor(sex_men), subset(NESARCw2_complex,black==1 | white==1), na.rm=T)
#svymean(~factor(black), subset(NESARCw2_complex,sex_men==1), na.rm=T)
#svymean(~factor(white), subset(NESARCw2_complex,sex_men==1), na.rm=T)

#hiv
CrossTable(NESARCw2_men$hiv_aids[NESARCw2_men$black==1 | NESARCw2_men$white==1])
svymean(~factor(hiv_aids), subset(NESARCw2_complex,black==1 | white==1), na.rm=T)


### ANALYSIS: DESCRIPTIVES by sex identity ###

describeBy(NESARCw2_men$age, NESARCw2_men$sex_identity_gay); t.test(NESARCw2_men$age ~ NESARCw2_men$sex_identity_gay)
CrossTable(NESARCw2_men$black, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$nationality, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$education, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$income, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$homeless_hx, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$kids, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$sex_identity, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$relationship, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$sex_yr, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$hiv_aids, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$std, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$abuse_sexual, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$abuse_physical, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$alcohol, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$alcohol_sex, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$drug_narcotic, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$drug_stimulant, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$drug_depressant, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NESARCw2_men$drug_hallucinogen, NESARCw2_men$sex_identity_gay, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)


### ANALYSIS: STRATIFIED by RACE ###

#partner gender

#all men
model = glm(hiv_aids~as.factor(sex_men)+as.factor(black)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=subset(NESARCw2_men,black==1 | white==1),family=binomial(link="logit"))
summary(model)
round(exp(coef(model)),2)
round(exp(confint(model)),2)

#black men
model = glm(hiv_aids~as.factor(sex_men)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=NESARCw2_black,family=binomial(link="logit"))
summary(model)
round(exp(coef(model)),2)
round(exp(confint(model)),2)

#white men
model = glm(hiv_aids~as.factor(sex_men)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=NESARCw2_white,family=binomial(link="logit"))
summary(model)
round(exp(coef(model)),2)
round(exp(confint(model)),2)

#test for interaction between MSM and black race
model1 = glm(hiv_aids~as.factor(sex_men)+as.factor(black)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=subset(NESARCw2_men,black==1 | white==1),family=binomial(link="logit"))
model2 = glm(hiv_aids~as.factor(sex_men)*as.factor(black)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=subset(NESARCw2_men,black==1 | white==1),family=binomial(link="logit"))
anova(model1,model2,test="Chisq")
rm(model1,model2)

#wald test for heterogeneity of effects; see PMID: 24022386
summary(glm(hiv_aids~as.factor(sex_men)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=NESARCw2_black,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(sex_men)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=NESARCw2_white,family=binomial(link="logit")))

#estimates and standard errors
b1 = 1.5060
b1_se = 0.5568
b2 = 2.7657
b2_se = 0.3884

#compute test statistic (Z-score)
zscore = (b2-b1)/(sqrt((b1_se^2) + (b2_se^2)))

#check against normal distribution
2*pnorm(-abs(zscore))


### COMPLEX SURVEY ESTIMATES ###

#black men

#some strata may have one sampled unit: http://r-survey.r-forge.r-project.org/survey/example-twostage.html
options(survey.lonely.psu="remove")

#specify survey design
NESARCw2_complex = svydesign(id=~psu, strata=~stratum, weights=~weight, data=NESARCw2)

#logistic regression using quasibinomial per: http://stackoverflow.com/questions/12953045/warning-non-integer-successes-in-a-binomial-glm-survey-packages
#partner
model = svyglm(hiv_aids~as.factor(sex_men)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug), family=quasibinomial, design=NESARCw2_complex)
summary(model)
round(exp(coef(model)),2)
round(exp(confint(model)),2)
AIC(model)

#non black men

#some strata may have one sampled unit: http://r-survey.r-forge.r-project.org/survey/example-twostage.html
options(survey.lonely.psu="remove")

#specify survey design
NESARCw2_complex = svydesign(id=~psu, strata=~stratum, weights=~weight, data=NESARCw2_white)

#logistic regression using quasibinomial per: http://stackoverflow.com/questions/12953045/warning-non-integer-successes-in-a-binomial-glm-survey-packages
#partner
model = svyglm(hiv_aids~as.factor(sex_men)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug), family=quasibinomial, design=NESARCw2_complex)
summary(model)
round(exp(coef(model)),2)
round(exp(confint(model)),2)
AIC(model)


### ANALYSIS: UNMEASURED CONFOUNDING ###

#naive black OR is 4.5; adjusted is 5.8
#naive white OR is 15.9 (log=2.8); adjusted is 11.5

#black men
model = glm(hiv_aids~as.factor(sex_men)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=NESARCw2_black,family=binomial(link="logit"))
summary(model)
round(exp(coef(model)),2)
round(exp(confint(model)),2)

#sensitivity analysis; see http://www.ncbi.nlm.nih.gov/pubmed/9750244
#g0 is the OR of U on HIV+
#p0 is prev of U in non-MSM
#p1 is prev of U in MSM

#naive
obsSensCCC(model, which=2, g0=seq(2,10,4), p0=seq(0,1,0.2), p1=seq(0,1,0.2), logOdds=F)

#informed from bayesian learning
obsSensCCC(model, which=2, g0=5, p0=seq(0,0.5,0.02), p1=seq(0,0.1,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=10, p0=seq(0,0.4,0.02), p1=seq(0,0.1,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=15, p0=seq(0,0.4,0.02), p1=seq(0,0.1,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=20, p0=seq(0,0.3,0.02), p1=seq(0,0.1,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=25, p0=seq(0,0.3,0.02), p1=seq(0,0.1,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=30, p0=seq(0,0.3,0.02), p1=seq(0,0.1,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=35, p0=seq(0,0.3,0.02), p1=seq(0,0.1,0.02), logOdds=F)

#white men
model = glm(hiv_aids~as.factor(sex_men)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=NESARCw2_white,family=binomial(link="logit"))
summary(model)
round(exp(coef(model)),2)
round(exp(confint(model)),2)

#sensitivity analysis; see http://www.ncbi.nlm.nih.gov/pubmed/9750244
#g0 is the OR of U on HIV+
#p0 is prev of U in non-MSM
#p1 is prev of U in MSM

#naive
obsSensCCC(model, which=2, g0=seq(2,10,4), p0=seq(0,1,0.2), p1=seq(0,1,0.2), logOdds=F)

#informed from bayesian learning
obsSensCCC(model, which=2, g0=2.5, p0=seq(0,0.1,0.02), p1=seq(0,0.8,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=5, p0=seq(0,0.1,0.02), p1=seq(0,0.4,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=10, p0=seq(0,0.1,0.02), p1=seq(0,0.3,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=15, p0=seq(0,0.1,0.02), p1=seq(0,0.3,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=20, p0=seq(0,0.1,0.02), p1=seq(0,0.3,0.02), logOdds=F)
obsSensCCC(model, which=2, g0=25, p0=seq(0,0.1,0.02), p1=seq(0,0.3,0.02), logOdds=F)


# ### PREVALENCE ###
# 
# #identity
# prev = sum(NESARCw2_black$sex_identity_gay, na.rm=T)/length(na.omit(NESARCw2_black$sex_identity_gay))
# prev_ci = 1.96 * sqrt((prev*(1-prev))/length(na.omit(NESARCw2_black$sex_identity_gay)))
# 
# #partner
# prev = sum(NESARCw2_black$sex_men, na.rm=T)/length(na.omit(NESARCw2_black$sex_men))
# prev_ci = 1.96 * sqrt((prev*(1-prev))/length(na.omit(NESARCw2_black$sex_men)))


### FIGURE ###

tiff("Figure1.tif",height=4,width=6,units='in',res=1200) 
pdf("Figure1.pdf",height=4,width=6) 

plot(x=c(0,2), y=c(0,18), type="l", xlim=c(0.2,0.9), ylim=c(0,10), axes=F, frame.plot=T, xlab="", ylab="POR(HIV | MSM)", lty=1)
lines(x=c(0,2), y=c(3.2,14.8), lty=5)
lines(x=c(0,2), y=c(-2.5,20.5), lty=4)
#plot(x=c(0,1), y=c(0,12), type="l", xlim=c(0.4,0.9), ylim=c(2,10), axes=F, frame.plot=T, xlab="", ylab="POR(HIV | MSM * Race)", lty=5)
#lines(x=c(0,1), y=c(0,6), lty=4)
#lines(x=c(0,1), y=c(0,9), lty=1)
axis(side=1, at=c(0.2,0.9), labels=c("non-MSM","MSM"))
#axis(side=2, at=c(0,2,4,6,8,10))
axis(side=4, at=c(0,2,4,6,8,10))
legend("topleft",lty=c(5,4,1),c("Black MSM ~ 5.8","White MSM ~ 11.5", "All MSM ~ 9.0"), cex=0.7)

dev.off() 