#################
# BAAMHS analysis
# Citation: Goldstein ND, Burstyn I, Welles SL. Bayesian approaches to racial disparities in HIV risk estimation among men who have sex with men. Epidemiology. 2017 Mar;28(2):215-220.
# 7/7/14 -- Neal Goldstein
#################


### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(car) #vif


### READ DATA ###

load("BAAMHS.RData")


### ANALYSIS: SN and SP ###

CrossTable(BAAMHS$sex_men[BAAMHS$hiv_aids==0], BAAMHS$any_anal[BAAMHS$hiv_aids==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(BAAMHS$sex_identity_gay[BAAMHS$hiv_aids==0], BAAMHS$any_anal[BAAMHS$hiv_aids==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(BAAMHS$sex_men[BAAMHS$hiv_aids==1], BAAMHS$any_anal[BAAMHS$hiv_aids==1], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(BAAMHS$sex_identity_gay[BAAMHS$hiv_aids==1], BAAMHS$any_anal[BAAMHS$hiv_aids==1], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


### ANALYSIS: ADJUSTED INFOMATIVE PRIOR ###

#crude
model = glm(hiv_aids~as.factor(any_anal),data=BAAMHS,family=binomial(link="logit"))
summary(model)
round(coef(model),2)
#variance is standard error (apprxs the standard deviation) squared
round(0.3892*0.3892,2)

#check for associations with hiv_aids, p<0.20
summary(glm(hiv_aids~age,data=BAAMHS,family=binomial(link="logit")))
#summary(glm(hiv_aids~as.factor(nationality),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(education),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(income),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(homeless_hx),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(kids),data=BAAMHS,family=binomial(link="logit")))
#summary(glm(hiv_aids~as.factor(religion),data=BAAMHS,family=binomial(link="logit")))
#summary(glm(hiv_aids~as.factor(relationship),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(std),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(abuse_sexual),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(abuse_physical),data=BAAMHS,family=binomial(link="logit")))
#summary(glm(hiv_aids~as.factor(alcohol),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(alcohol_sex),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(recent_drug),data=BAAMHS,family=binomial(link="logit")))

#check for associations with any_msm_behavior, p<0.10 and p<0.10 with hiv_aids
#summary(glm(any_anal~age,data=BAAMHS,family=binomial(link="logit")))
#summary(glm(any_anal~as.factor(nationality),data=BAAMHS,family=binomial(link="logit")))
#summary(glm(any_anal~as.factor(education),data=BAAMHS,family=binomial(link="logit")))
summary(glm(any_anal~as.factor(income),data=BAAMHS,family=binomial(link="logit")))
summary(glm(any_anal~as.factor(homeless_hx),data=BAAMHS,family=binomial(link="logit")))
#summary(glm(any_anal~as.factor(kids),data=BAAMHS,family=binomial(link="logit")))
#summary(glm(any_anal~as.factor(religion),data=BAAMHS,family=binomial(link="logit")))
summary(glm(any_anal~as.factor(relationship),data=BAAMHS,family=binomial(link="logit")))
summary(glm(any_anal~as.factor(std),data=BAAMHS,family=binomial(link="logit")))
summary(glm(any_anal~as.factor(abuse_sexual),data=BAAMHS,family=binomial(link="logit")))
summary(glm(any_anal~as.factor(abuse_physical),data=BAAMHS,family=binomial(link="logit")))
summary(glm(any_anal~as.factor(alcohol),data=BAAMHS,family=binomial(link="logit")))
summary(glm(any_anal~as.factor(alcohol_sex),data=BAAMHS,family=binomial(link="logit")))
summary(glm(any_anal~as.factor(recent_drug),data=BAAMHS,family=binomial(link="logit")))

#full model with potential confounders
summary(glm(hiv_aids~as.factor(any_anal)+as.factor(income)+as.factor(homeless_hx)+as.factor(std)+as.factor(abuse_sexual)+as.factor(abuse_physical)+as.factor(alcohol_sex)+as.factor(recent_drug),data=BAAMHS,family=binomial(link="logit")))

#check for multicollinearity, VIF>=10
vif(glm(hiv_aids~as.factor(any_msm_behavior)+age+as.factor(income)+as.factor(homeless_hx)+as.factor(std)+as.factor(abuse_sexual)+as.factor(abuse_physical)+as.factor(alcohol_sex)+as.factor(drug_narcotic)+as.factor(drug_stimulant)+as.factor(drug_hallucinogen),data=BAAMHS,family=binomial(link="logit")))

#backward remove nonsignificant vars and check for change in estimate >10%
summary(glm(hiv_aids~as.factor(any_anal)+as.factor(income)+as.factor(homeless_hx)+as.factor(std)+as.factor(abuse_sexual)+as.factor(abuse_physical)+as.factor(alcohol_sex)+as.factor(recent_drug),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(any_anal)+as.factor(homeless_hx)+as.factor(std)+as.factor(abuse_sexual)+as.factor(abuse_physical)+as.factor(alcohol_sex)+as.factor(recent_drug),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(any_anal)+as.factor(homeless_hx)+as.factor(std)+as.factor(abuse_sexual)+as.factor(alcohol_sex)+as.factor(recent_drug),data=BAAMHS,family=binomial(link="logit")))
summary(glm(hiv_aids~as.factor(any_anal)+as.factor(homeless_hx)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=BAAMHS,family=binomial(link="logit")))

#final model w/ same covariates from NESARC
model = glm(hiv_aids~as.factor(any_anal)+as.factor(relationship)+as.factor(std)+as.factor(abuse_sexual)+as.factor(recent_drug),data=BAAMHS,family=binomial(link="logit"))
summary(model)
round(coef(model),2)
#variance is standard error (apprxs the standard deviation) squared
round(0.44682*0.44682,2)

round(exp(coef(model)),2)
round(exp(confint(model)),2)


### PREVALENCE ###

#gold standard (anal sex)
prev = sum(BAAMHS$any_anal, na.rm=T)/length(na.omit(BAAMHS$any_anal))
prev_ci = 1.96 * sqrt((prev*(1-prev))/length(na.omit(BAAMHS$any_anal)))

#identity
prev = sum(BAAMHS$sex_identity_gay, na.rm=T)/length(na.omit(BAAMHS$sex_identity_gay))
prev_ci = 1.96 * sqrt((prev*(1-prev))/length(na.omit(BAAMHS$sex_identity_gay)))

#partner
prev = sum(BAAMHS$sex_men, na.rm=T)/length(na.omit(BAAMHS$sex_men))
prev_ci = 1.96 * sqrt((prev*(1-prev))/length(na.omit(BAAMHS$sex_men)))

