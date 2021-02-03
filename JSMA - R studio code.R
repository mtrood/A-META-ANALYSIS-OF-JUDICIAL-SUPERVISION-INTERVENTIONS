#############################################################################################
# Judicial Supervision Meta-analysis: Analysis
#############################################################################################
########################################
### I) Table of Contents 
########################################
#   I:  Table of Contents
#  II:  References
#   1:  Packages
#   2:  Data cleaning and wrangling
#       2.1: Data
#       2.2: Convert relevant variables to numeric
#       2.3: Incidents subset
#       2.4: Descriptives/data checking
#   3:  Calculate the Relative Incident Rate Ratio 
#   4:  Weighted Mean Log Relative Incident Rate Ratio (RIRR): Fixed Effects Model
#   5:  Random effects model: DerSimonian-Laird Estimate 
#   6:  Multiplicative Variance Adjustment
#       6.1:  Multiplicative Variance Adjustment: Fixed Effects Model
#       6.2:  Sensitivity analysis: Multiplicative Variance Adjustment: Fixed Effects Model
#             6.2.1: Alternative over-dispersion factors
#	            6.2.2: Publication Bias
#   7:  Forest Plot
#   8:  Fixed effect model with categorical (Multiplicative Variance Adjustment) and continuous moderators
#       8.1: Methodological 
#             8.1.1:      Year of publication
#             8.1.2:      Retrospective vs Contemporary studygroups
#             8.1.3:      Pre and post period (months)
#             8.1.4:      Months measured post court-program
#             8.1.5:      Published/unpublished
#             8.1.6:      Published by NPC Research 
#       8.2: Bias and quality 
#             8.2.1:      Score on the modified Maryland Scientific Methods Scale
#             8.2.2:      ROBINS confounding  
#             8.2.3:      ROBINS selection of participants rating - Moderate and serious
#             8.2.4:      ROBINS classification of treatment status rating - Low and moderate
#             8.2.5:      ROBINS deviation from intended intervention
#             8.2.6.1:    ROBINS missing rating - Low and moderate
#             8.2.6.2:    ROBINS missing rating - Low and serious
#             8.2.7:      ROBINS measurement rating - Low and moderate
#             8.2.8:      ROBINS selection of reported result 
#       8.3: Sample characteristics  
#             8.3.1:      Lifetime arrests - Treatment
#             8.3.2:      Lifetime arrests - Control
#             8.3.3:      Mean age - Treatment
#             8.3.4:      Mean age - Control
#             8.3.5:      Percent male - Treatment
#             8.3.6:      Percent male - Control
#             8.3.7:      Percent white - Treatment 
#             8.3.8:      Percent white - Control 
#             8.3.9:      Percent program graduates - Treatment
#             8.3.10:     Percent program terminates - Treatment
#             8.3.11:     Percent treatment completed - Treatment
#             8.3.12:     Percent treatment completed - Control
#       8.4: Programmatic 
#             8.4.1.1:    Court type - Drug court
#             8.4.1.2:    Court type - Driving-under-the-influence court
#             8.4.1.3:    Court type - Family drug court
#             8.4.2:      Minimum program duration (months) 
#             8.4.3:      Average program duration (months)
#             8.4.4:      Accepted violent offenders  
#             8.4.5:      Reported accepting violent offenders  
#             8.4.6:      Accepted felons
#             8.4.7:      Reported accepting felons
#             8.4.8:      Reported that charges were expunged
#             8.4.9:      Stage of disposition
#             8.4.10:     Single vs multiple judges presiding over court hearings
#             8.4.11:     Reported that a single judge presided over court hearings
#             8.4.12:     Reported team training
#             8.4.13:     Reported training judge
#             8.4.14:     Reported training judge or team
#             8.4.15:     Average minutes per court review hearing reported for court
#             8.4.16:     Reported frequent team meetings
#             8.4.17:     Reported using sanctions 
#             8.4.18:     Reported using graduated sanctions
#             8.4.19.1:   Reported sanctions or graduated sanctions - Graduated 
#             8.4.19.2:   Reported sanctions or graduated sanctions - Sanctions 
#             8.4.20:     Reported using rewards 
#             8.4.21:     Reported using sanctions and rewards 
#             8.4.22.1:   Frequency of review hearings - Weekly
#             8.4.22.2:   Frequency of review hearings - Matched to level of risk and need
#             8.4.23:     Frequent review hearings collapsed
#       8.5: Treatment
#             8.5.1:      Reported measuring risk - Treatment
#             8.5.2:      Reported measuring risk - Control
#             8.5.3:      Reported individualising treatment - Treatment
#             8.5.4:      Reported individualising treatment - Control
#             8.5.5:      Reported drug or alcohol treatment - Treatment
#             8.5.6:      Reported drug or alcohol treatment - Control
#             8.5.7:      Reported drug or alcohol testing - Treatment
#             8.5.8:      Reported drug or alcohol testing - Control
#             8.5.9:      Reported psychological treatment - Treatment
#             8.5.10:     Reported psychological treatment - Control
#             8.5.11:     Reported counselling - Treatment
#             8.5.12:     Reported counselling - Control
#             8.5.13:     Reported support groups - Treatment
#             8.5.14:     Reported support groups - Control
#             8.5.15:     Reported probation - Control
#             8.5.16:     Reported other treatment - Aftercare - Treatment
#             8.5.17:     Reported other treatment - Relapse prevention - Treatment
#             8.5.18:     Reported other treatment - Self-help - Treatment
#             8.5.19:     Total treatments required - Treatment
#             8.5.20:     Total treatments required - Control
#             8.5.21:     Total treatments offered - Treatment
#             8.5.22:     Total treatments offered - Control
#   9:  Z test for subgroup differences (categorical moderators) 
#       9.1: Multiplicative Variance Adjustment  
#       9.2: Z test 
#  10:  Benjamini-Hochberg  adjustments 
#       10.1: Categorical moderators 
#       10.2: Continuous moderators
#  11:  Pooled effects table
#  12:  Categorical moderator table
#       12.1: Format table 
#       12.2: Add text, move and merge columns, and export
#  13:  Continuous moderator table
#       13.1: Amend figures
#       13.2: Edit text and export
#  14:  Health and well-being outcomes
#       14.1: Separate target rows 
#       14.2: Calculate effect sizes and statistical significance 
########################################
### II) References
########################################

### All mathematical formula were drawn from the following sources:
### 
### Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2011). Introduction to meta-analysis. John Wiley & Sons.
###     
### Deeks, J. J., Higgins, J. P., Altman, D. G., & Cochrane Statistical Methods Group. (2008). Analysing data and undertaking  
###     meta‐analyses. In J. P. Higgins & D. G. Altman (Eds.), Cochrane handbook for systematic reviews of interventions (pp. 243-
###     296). Chichester, UK: John Wiley & Sons, Ltd.
###     
### DerSimonian, R., & Kacker, R. (2007). Random-effects model for meta-analysis of clinical trials: an update. Contemporary 
###     clinical trials, 28(2), 105-114.
###
### Duval, S., & Tweedie, R. (2000). A Nonparametric "Trim and Fill" Method of Accounting for Publication Bias in Meta-Analysis. 
###     Journal of the American Statistical Association, 95(449), 89-98. doi:10.1080/01621459.2000.10473905
###     
### Higgins, J. P. and Thompson, S. G. (2002). "Quantifying heterogeneity in a meta-analysis." Statistics in medicine 21(11): 
###     1539-1558.
###
### Higgins, J. P., Thompson, S. G., Deeks, J. J., & Altman, D. G. (2003). Measuring inconsistency in
###     meta-analyses. Bmj, 327(7414), 557-560.
###
### Jones, HE. 2005. Measuring Effect Size in Area-Based Crime Prevention Research.
###     Unpublished M.Phil. thesis. Cambridge, UK: Statistical Laboratory, Cambridge University
###
### Viechtbauer, W. (2010). "Metafor: meta-analysis package for R." R package version 2010: 1-0.
###
### Welsh, B. C., & Farrington, D. P. (2007). Closed-circuit television surveillance. In B. C. Welsh & D. P. Farrington (Eds.), 
###     Preventing crime: What works for children, offenders, victims and places (pp. 193-208). Springer New York. 
###     https://doi.org/10.1007/978-0-387-69169-5_13 
###
### Welsh, B. C., & Farrington, D. P. (2008). Effects of closed circuit television surveillance on crime.
###     Campbell systematic reviews, 4(1), 1-73.
###
########################################
### 1) Packages 
########################################
library(plyr)
library(irr)
library(dplyr)
library(tidyverse)
library(stats)
library(metafor)
library(R.utils)
library(forestplot)
library(lattice)
library(lubridate)
library(stringr)
library(data.table)
library(sjstats)
library(ggplot2)
library(stringi)
########################################
### 2) Data cleaning and wrangling
########################################
########################################
### 2.1) Data
########################################
# Import data
MA_EX <- read.csv(url("https://raw.githubusercontent.com/mtrood/A-META-ANALYSIS-OF-JUDICIAL-SUPERVISION-INTERVENTIONS/main/JSMA.csv"))

# Replace NRs with NA

MA_EX[ MA_EX == "NR" ] <- NA
########################################
### 2.2) Convert relevant variables to numeric
########################################
MA_EX_PT <- MA_EX
# Convert appropriate factor variables to numeric 
MA_EX_PT_num<-MA_EX_PT%>%select(starts_with("LT"), starts_with("mean"),starts_with("pt"), -"pt_begins",
                                starts_with("male"),starts_with("white"),starts_with("bl"),
                                contains("program"),contains("Average"), starts_with("total"))%>%
                                apply(2,function(x) as.numeric(as.character(x))) 

MA_EX_PT_no_num<-MA_EX_PT%>%select(-starts_with("LT"), -starts_with("mean"),-starts_with("pt"), "pt_begins",
                                   -starts_with("male"),-starts_with("white"),-starts_with("bl"),
                                   -contains("program"),-contains("Average"), -starts_with("total")) #subset dataframe of factor variables

MA_EX_PT_final<-cbind(MA_EX_PT_num,MA_EX_PT_no_num) #bind numeric and factor variables
MA_EX_PT_final<-MA_EX_PT_final[,names(MA_EX_PT)]#reset order of variables

##Round Incidents to nearest whole number

MA_EX_PT_final$bl_in_rounded_treat <- round(MA_EX_PT_final$bl_incidents_treat)
MA_EX_PT_final$bl_in_rounded_cont <-  round(MA_EX_PT_final$bl_incidents_cont)
MA_EX_PT_final$pt_in_rounded_treat <- round(MA_EX_PT_final$pt_incidents_treat)
MA_EX_PT_final$pt_in_rounded_treat <- ifelse(MA_EX_PT_final$pt_in_rounded_treat==0, 1, MA_EX_PT_final$pt_in_rounded_treat)
MA_EX_PT_final$pt_in_rounded_cont <-  round(MA_EX_PT_final$pt_incidents_cont)

########################################
### 2.3) Incidents subset 
########################################
###Separate only outcomes with "serious" rating on ROBINS
incidents <-MA_EX_PT_final%>%filter(overall_robins=="serious")
# incidents$Measure <-droplevels(incidents$Measure) 

#Separate only those with an outcomes that have arrest counts
incidents_bound <-incidents%>%filter(Measure=="mean arrests")


#Exclude studies because they contain duplicate effects to other studies 
incidents_bound<-incidents_bound%>%filter(!(studyid==11|      # Same sample as Study ID 12 with shorter pre/post measurement periods
                                              studyid==33|    # Identical sample and effect size to Study ID 32
                                              studyid==50))   # Same sample as STudy ID 57 but with missing data.

#To re-order all rows ascending from study ID 
incidents_bound_final <- incidents_bound %>% arrange(studyid)

########################################
### 2.4) Descriptives/data checking
########################################
percent_check<-incidents_bound%>% 
  apply(.,2,function(x) as.character(x))%>%
  apply(.,2,function(x) str_detect(x,"%")) #check for percent (%) symbols and remove

#how many % symbols observed
sum(percent_check==TRUE, na.rm=T)

########################################
### 3) Calculate the Relative Incident Rate Ratio (RIRR)
########################################
### Jones (2005, p.7) calculates the RIRR using a poisson distribution using the following, where r is the event-count 
###                   and n is the sample number.     
###               
###               Pre         Post          
### Experimental: r11/n11     r12/n12
### Comparison  : r21/n21     r22/ n22

# The relative incident rate ratio (RIRR) is calculated as follows:
#                          ((r.11/n.11)*(r.22/n.22))/((r.12/n.12)*(r.21/n.21))
# which is equivalent to:  ((baseline_treatment_mean)*(post-test_control_mean)/
#                          (post-test_treatment_mean)*(baseline_control_mean))

incidents_bound_final$RIRR <- ((incidents_bound_final$bl_treat_m)*(incidents_bound_final$pt_cont_m))/
                              ((incidents_bound_final$pt_treat_m)*(incidents_bound_final$bl_cont_m))
###Reorder the data-frame in terms of the size of the RES 
incidents_bound_final <- incidents_bound_final[order(incidents_bound_final$RIRR),]

RIRR <- incidents_bound_final$RIRR

#Variance: Jones (2005, p16-17) calculates the variance of natural log relative incident rate ratio using the incidents 
#          counts as follows: 
#          VAR(LRIRR) = (1/a) + (1/b) + (1/c) + (1/d). 
#          where a is the incident count for the treatment group at time 1, b is the incident count for the treatment group at time 2,
#          where c is the incident count for the control group at time 1, d is the incident count for the control group at time 2,

incidents_bound_final$RIRR_var <- (1 / incidents_bound_final$bl_in_rounded_treat) + 
                                 (1 / incidents_bound_final$pt_in_rounded_treat) +
                                 (1 / incidents_bound_final$bl_in_rounded_cont) + 
                                 (1 / incidents_bound_final$pt_in_rounded_cont)  

v <- incidents_bound_final$RIRR_var

# Note: Farrington & Welsh (2013, p.5) note that several authors have multiplied effect size's variance to account 
# for over-dispersion, usually by doubling the variance of each effect size (e.g., Welsh & Farrington, 2007). 
# This can be achieved by multiplying 'v' by 2. 
v <- v*2

incidents_bound_final$RIRR_var <- incidents_bound_final$RIRR_var *2
### Log(RIRR)
LRIRR <- log(RIRR)
### Z score: RIRR
z<-LRIRR/sqrt(v)
### Z score: p value
p<- round(2 * pnorm(abs(z), lower.tail = FALSE), digits = 3)

### Z score confidence interval 
lower.CI<-exp(LRIRR-1.96*sqrt(v))
upper.CI<-exp(LRIRR+1.96*sqrt(v))


### Create a histogram for each study's RES. 
histogram_df <- as.data.frame(cbind(RIRR, v))
es_histogram <- ggplot(histogram_df, aes(x=RIRR)) + 
                geom_histogram(aes(y=..density..), color="darkblue", fill="lightblue", binwidth= 0.1)+ 
                geom_density(alpha=.2, fill="#FF6666") 

########################################
### 4) Weighted Mean Log Relative Incident Rate Ratio (RIRR): Fixed Effects Model
########################################
### populationLRIRR<-(sum(LRIRR/v))/(sum(1/v))
fe_poplnLRIRR<-(sum(LRIRR/v))/(sum(1/v))
fe_poplnLRIRR

### Exponential of fe_poplnLRIRR
fe_poplnRIRR<-exp(fe_poplnLRIRR)
fe_poplnRIRR

### Variance poplnRIRR
fe_var_poplnLRIRR<-1/sum(1/v)
fe_var_poplnLRIRR

### Standard Error poplnRIRR 
### According to (Borenstein et al., 2011, p. 37), an estimate of the population standard error is given 
### by SElogOddsratio = sqrt(LogOddsRatioVariance). We will take the variance on the RIRR using the same calculation. 
fe_se_poplnLRIRR<- sqrt(fe_var_poplnLRIRR)
fe_se_poplnLRIRR  

### Calculate the 'test statistic' cited by Jones (2005)
fe_teststat <- fe_poplnLRIRR/sqrt(fe_var_poplnLRIRR)
fe_teststat
### Test statistic P value
fe_pvalue <-2 * pnorm(abs(fe_teststat), lower.tail = FALSE)
fe_pvalue 

### Population RIRR confidence interval
fe_RIRR_CI<-c(exp(fe_poplnLRIRR-1.96*sqrt(fe_var_poplnLRIRR)),
                  exp(fe_poplnLRIRR+1.96*sqrt(fe_var_poplnLRIRR)))
fe_RIRR_CI

### Q statistic for homogeneity (Jones, 2005)
Q <- sum((LRIRR-fe_poplnLRIRR)^2/v)
Q

#Borenstein et al. (2011, p.12) note that Q follows a chi-squared distribution, and that a p-value can be
#calculated on Q using k-1 degrees of freedmon
Q_p<- pchisq(Q, df=(length(LRIRR)-1), lower.tail = FALSE)
Q_p

###I^2 - Higgins et al. (2003) calculate I^2 as follows:
###      ((Q-df)/Q)*100 
fe_I <- ((Q-(length(LRIRR)-1))/Q)*100 
fe_I

###Create a vector with the row name, Relative Incident Rate Ratio, Confidence Interval, Z statistic, Z p-value, 
###Credibility Interval, Q statistic and I2
pooled_effect_fe<-c('fixed effect model', fe_poplnRIRR, fe_RIRR_CI, fe_teststat, 
                    fe_pvalue, '-','-', Q, Q_p, fe_I)

###Add the names for each column
names(pooled_effect_fe)<-c("Model", "Relative Incident Rate Ratio", "Conlow", "Conup", "Test statistic", "p-value",
                           "credlow", "credup", "Q statistic","Q_pv", "I2")
pooled_effect_output<- pooled_effect_fe

########################################
### 5) Random effects model: DerSimonian-Laird Estimate
########################################
### We will follow Borenstein et al.'s (2011, pp. 72-3) equations and definitions for calculating Tau^2 (T.sq hereafter) and ultimately the DerSimonian-Laird Estimate. 
### T.sq is estimated as follows (equation 12.2: T.sq = (Q-df)/C 
###   Where C = the sum of the study weights minus the sum of each study weight squared and divided by 
###   sum of the the study weights, or C = sum(1/v)-(sum((1/v)^2)/sum(1/v))
dle_C <- sum(1/v)-(sum((1/v)^2)/sum(1/v))

### calculate T.sq
dle_T.sq<- (Q-(length(LRIRR)-1))/dle_C
dle_T.sq

### DerSimonian-Laird study weights 1/ V*yi 
### Where V*yi = Vyi +T.sq
dle_v <- v + dle_T.sq
dle_v
### DerSimonian-Laird LRIRR
dle_poplnLRIRR<- (sum(LRIRR/dle_v ))/(sum(1/dle_v))
dle_poplnLRIRR

### DerSimonian-Laird LRIRR exponentiated 
dle_poplnRIRR<-exp(dle_poplnLRIRR)
dle_poplnRIRR

### DerSimonian-Laird LRIRR Variance
dle_var_poplnLRIRR<-1/sum(1/dle_v)
dle_var_poplnLRIRR

### DerSimonian-Laird test statistic 
dle_teststat<-dle_poplnLRIRR/sqrt(dle_var_poplnLRIRR)
dle_teststat

### DerSimonian-Laird test statistic p-value
dle_pvalue<-2 * pnorm(abs(dle_teststat), lower.tail = FALSE)
dle_pvalue

### DerSimonian-Laird test statistic CI
dle_teststat_CI<-c(exp(dle_poplnLRIRR-1.96*sqrt(dle_var_poplnLRIRR)),
          exp(dle_poplnLRIRR + 1.96*sqrt(dle_var_poplnLRIRR)))

### Calculate Q
### Note: we will calculate Q using code from Viectbauer's (2010) 'rma.uni' function's source code. This will produce the same Q statistic as the fixed-effect model. 
dle_Y <- as.matrix(LRIRR)
dle_W <- diag(1/v)
dle_X <- as.matrix(rep(1, length(LRIRR)))
dle_P <- dle_W - dle_W %*% dle_X %*% solve(
  t(dle_X) %*% dle_W %*% dle_X) %*% t(
    dle_X) %*% dle_W 

dle_Q <- max(0, c(crossprod(dle_Y, dle_P) %*% dle_Y))

###Calculate a p-value for Q
dle_Q_p<- pchisq(dle_Q, df=(length(LRIRR)-1), lower.tail = FALSE)
dle_Q_p

### T (our estimate of the true standard deviation) is simply the square root of T.sq (Borenstein, et al, 2011, p.116)
dle_T<-sqrt(dle_T.sq)
dle_T

### Calculate I^2
### We will follow Viectbauer's (2010) approach to calculating I^2 using Higgins and Thompson's (2002, p.1546) formula (9)
### I^2 = Tau^2 / (Tau^2 + Vsig)
###                 where Vsig =  Sum(wi)*k-1/ sum(wi)^2 - sum(wi^2)
Vsig <- sum(1/v*(length(LRIRR)-1)) / (sum(1/v)^2 - sum((1/v)^2))
dle_I2 <- dle_T.sq/(dle_T.sq + Vsig) * 100
dle_I2

###Credibility Interval 
###Borenstein et al. (2011, p. 129) cite Higgins et al.'s (year of study not reported) formulae for
###calculating the credibility interval
###   LLpred = M* - tdf * (sqrt(T.sq + Vm*))
###   ULpred = M* + tdf * (sqrt(T.sq + Vm*))
###     where M* = mean effect size in the sample, tdf is the t-value corresponding to the desired alpha level
###     with the relevant degrees of freedom, and Vm* is the variance of M*

#Calculate tdf. 
# qt(1-alpha/2, k-2)
dle_tdf <- qt(1-0.05/2, (length(LRIRR)-2))
dle_tdf

dle_LLcred <- dle_poplnLRIRR - dle_tdf*(sqrt(dle_T.sq + dle_var_poplnLRIRR))
dle_ULcred <- dle_poplnLRIRR + dle_tdf*(sqrt(dle_T.sq + dle_var_poplnLRIRR))
dle_CI <- c(exp(dle_LLcred), exp(dle_ULcred))
dle_CI

###Create a vector with the row name, Relative Incident Rate Ratio, Confidence Interval, Z statistic, Z p-value, 
###Credibility Interval, Q statistic and I2
pooled_effect_dle<-c('DerSimonian-Laird Estimate: Random effects model', dle_poplnRIRR, dle_teststat_CI, 
                     dle_teststat, dle_pvalue, dle_CI, dle_Q, dle_Q_p, dle_I2)

###Add the names for each column
names(pooled_effect_dle)<-c("Model", "Relative Incident Rate Ratio", "Conlow", "Conup", "Test statistic", "p-value",
                            "credlow", "credup", "Q statistic","Q_pv", "I2")
pooled_effect_output<- bind_rows(pooled_effect_output, pooled_effect_dle)
########################################
### 6) Multiplicative Variance Adjustment
########################################
########################################
### 6.1) Multiplicative Variance Adjustment: Fixed Effects Model 
########################################
### Note: This section employs Jones'(2005) formulae and code and Farrington's (2013) Formula for Vm
###MVA Population LogRIRR
mva_PLRIRR <- fe_poplnLRIRR
mva_PLRIRR

### Note: that the estimate is in log units. Exponentiate the RIRR
mva_RIRR <- exp(fe_poplnLRIRR)
mva_RIRR

###Note that we calculate the Q and I^2 (I2) statistics using the same inputs as the fixed effects model.
###Importantly, this means using VAR(LRIRR) to calculate Q and not the adjusted variance (va). 
###Q statistic 
Q_mva <-sum((LRIRR-mva_PLRIRR)^2/v)
Q_mva

#Borenstein et al. (2011, p.12) note that Q follows a chi-squared distribution, and that a p-value can be
#calculated on Q using k-1 degrees of freedom
Q_mva_p<- pchisq(Q_mva, df=(length(LRIRR)-1), lower.tail = FALSE)
Q_mva_p

### Adjusted Variance
### Farrington & Welsh (2013, p.9) describe the following formula for the multiplicative variance adjustment
### Vm = Vf * Ve
### Where Vm = variance in MVA model
###       Vf = variance in the fixed effects model
###       Ve = extra variance where Ve = Q / df 

Ve <- (Q_mva/(length(LRIRR)-1)) 
Vm <- round(fe_var_poplnLRIRR*Ve, digits = 10) 
Vm 
### Jones (2005, p28) calculated Vm by squaring the standard error of the intercept estimated by a 
### linear model. To test if they are the same value:
### Linear Model 
model_lm <-lm(LRIRR ~ 1,weights = 1/v)
summary(model_lm)

Vm_alt<- round(summary(model_lm)$coefficients[1, 2]^2, digits = 10)
Vm_alt==Vm
### For the sake of brevity, Jones' (2005, p28) formula will be used in subsequent sections

### Note: Jones (2005, p28) calculates the adjusted 95% CI for the population odds-ratio as follows
### using the standard error of the intercept (see 'Adjusted Variance' section below for further details). 
mva_RIRR_CI_low <- exp(mva_PLRIRR-1.96*sqrt(Vm))
mva_RIRR_CI_upp <- exp(mva_PLRIRR+1.96*sqrt(Vm))
mva_RIRR_CI <-  c(exp(mva_PLRIRR-1.96*sqrt(Vm)), exp(mva_PLRIRR+1.96*sqrt(Vm)))
mva_RIRR_CI

### Calculate Jones' (2005) test-statistic. 
mva_ts <-mva_PLRIRR/sqrt(Vm)
mva_ts
### P value
mva_ts_p <-2 * pnorm(abs(mva_ts), lower.tail = FALSE)
mva_ts_p

###I^2 - Higgins et al. (2003) calculate I^2 as follows:
###I^2 = ((Q-df)/Q)*100
I_mva <- ((Q_mva-(length(LRIRR)-1))/Q_mva)*100  
I_mva

###Create a vector with the row name, Relative Incident Rate Ratio, Confidence Interval, Z statistic, Z p-value, 
###Credibility Interval, Q statistic and I2
pooled_effect_mva<-c('Multiplicative Variance Adjustment', mva_RIRR, mva_RIRR_CI_low,mva_RIRR_CI_upp,
                     mva_ts, mva_ts_p, '-','-', Q_mva, Q_mva_p, I_mva)

###Add the names for each column
names(pooled_effect_mva)<-c("Model", "Relative Incident Rate Ratio", "Conlow", "Conup", "Test statistic", "p-value",
                            "credlow", "credup", "Q statistic","Q_pv", "I2")

pooled_effect_output<- bind_rows(pooled_effect_output, pooled_effect_mva)

########################################
### 6.2) Sensitivity analysis: Multiplicative Variance Adjustment
########################################
########################################
### 6.2.1) Alternative over-dispersion factors
########################################
### Code a dataframe that includes the model statistics for every possible combination of over-disperion factor from 1 to 10 by 0.01
iterative_function <- function(vector_yi, vector_vi){
  overdispersion_factor <-  (seq(from=1, to=10, by=0.01)) # create a vector with the dispersion factors values
  vi_long <- as.data.frame(do.call(cbind, replicate(901, vector_vi, simplify = FALSE)))
  vi_long_adjusted <- t(t(vi_long)*overdispersion_factor) # multiply each vi by column by the successive overdispersion factor 
  column_combinations <-cbind(seq(from=1, to=901)) # create a df that will allow us to index the vi columns when we reiterate the model
  models_list <- apply(column_combinations, 1, function(w){
    lm(log(vector_yi) ~ 1 ,weights = 1/vi_long_adjusted[, w[1]]) # gives us a list of the models
  })
  summary_data <- lapply(models_list, function(x)summary(x)$coefficients[1,1:2]) # extracts each model's coefficient and E
  summary_df <- data.frame(matrix(unlist(summary_data), nrow=length(summary_data), byrow=T)) # converts the list into a df
  colnames(summary_df) <- c("PLRIRR", "PLRIRR_se") 
  summary_df$Q <- apply(column_combinations, 1, function(y){ 
    sum((log(vector_yi)- summary_df[y[1], 1])^2/vi_long_adjusted[, y[1]])                                   
  })
  summary_df$Q_p <- pchisq(summary_df$Q, df=(length(vector_yi)-1), lower.tail = FALSE)
  summary_df$Ve <- (summary_df$Q/(length(vector_yi)-1)) 
  summary_df$Vm <- apply(column_combinations, 1, function(z){
    1/sum(1/vi_long_adjusted[, z[1]])*summary_df[z[1], 5]
  })
  summary_df$Z <- summary_df$PLRIRR/sqrt(summary_df$Vm)
  summary_df$Z_p <- 2 * pnorm(abs(summary_df$Z), lower.tail = FALSE)
  summary_df$ci.lb <- summary_df$PLRIRR-1.96*sqrt(summary_df$Vm)
  summary_df$ci.ub <- summary_df$PLRIRR+1.96*sqrt(summary_df$Vm)
  summary_df$I <- ((summary_df$Q-(length(vector_yi)-1))/summary_df$Q)*100  
  summary_df[c(1,9:10)] <- lapply(summary_df[c(1,9:10)], exp)
  summary_df$overdispersion_factor <- overdispersion_factor
  return(summary_df)
}

unadjusted_RIRR_var <- incidents_bound_final$RIRR_var/2 # Divide the variance by two because it was doubled in Section 3

sensitivity_analysis <- iterative_function(incidents_bound_final$RIRR, unadjusted_RIRR_var)

### Call the lowest non-significant p-value for Q 
sensitivity_analysis[which.max(sensitivity_analysis$Q_p >=0.05), ]

### Call the  highest I2 equal to or below 0 
sensitivity_analysis[which.max(sensitivity_analysis$I <=0.001), ]

########################################
### 6.2.1) Publication Bias
########################################
### Note: we will employ Duval and Tweedie's (2000) trim and fill procedure on a fixed effect model using Viectbauer's (2010)
###       trimfill.rma.uni function and will perform the multiplicative variance adjustment on the resulting fixed effect model.
###       This is possible given the procedure doesn't incorporate the population standard error or population variance
###       until the final model in calculated on the trimmed and filled data (Duval and Tweedie's, 2000, p. 94).

### First create an rma.uni object using the LRIRR and LRIRR variance values in the yi and vi argruments
trim_and_fill_incidents_bound_final <- incidents_bound_final
trim_and_fill_incidents_bound_final$LRIRR <- log(trim_and_fill_incidents_bound_final$RIRR)

trim_and_fill_rma <- rma(trim_and_fill_incidents_bound_final$LRIRR, trim_and_fill_incidents_bound_final$RIRR_var, 
                                          method = 'FE', data = trim_and_fill_incidents_bound_final)

### Run metafor's trimfill on the rma object
trim_and_fill_model <- trimfill(trim_and_fill_rma)

### Funnel plot including the estimated missing values
dev.new(width=22, height=22)

funnel(trim_and_fill_model, legend = TRUE, atransf=exp, xlab = 'RIRR')

# Trim-and-fill population LRIRR and RIRR
trim_and_fill_LRIRR <- trim_and_fill_model$beta[1]
trim_and_fill_LRIRR

trim_and_fill_RIRR <- round(exp(trim_and_fill_LRIRR), digits = 2)
trim_and_fill_RIRR

### Perform the MVA
trim_and_fill_Ve <- trim_and_fill_model$QE/(length(trim_and_fill_model$yi)-1)
trim_and_fill_Vm <- round(trim_and_fill_model$se^2*trim_and_fill_Ve, digits = 10) 
trim_and_fill_Vm

### Calculate the adjusted confidence interval 
trim_and_fill_RIRR_CI_low <- exp(trim_and_fill_LRIRR-1.96*sqrt(trim_and_fill_Vm))
trim_and_fill_RIRR_CI_upp <- exp(trim_and_fill_LRIRR+1.96*sqrt(trim_and_fill_Vm))
trim_and_fill_RIRR_CI <-  c(exp(trim_and_fill_LRIRR-1.96*sqrt(trim_and_fill_Vm)), exp(trim_and_fill_LRIRR+1.96*sqrt(trim_and_fill_Vm)))
trim_and_fill_RIRR_CI

### Calculate Jones' (2005) test-statistic. 
trim_and_fill_ts <-trim_and_fill_LRIRR/sqrt(trim_and_fill_Vm)
trim_and_fill_ts
### P value
trim_and_fill_ts_p <-2 * pnorm(abs(trim_and_fill_ts), lower.tail = FALSE)
round(trim_and_fill_ts_p, digits = 5)

### Rosenthal's Fail Safe N
rosenthal_fsn<- fsn(yi = log(incidents_bound_final$RIRR), vi = incidents_bound_final$RIRR_var)
rosenthal_fsn

### Orwin's Fail Safe N 
orwin_fsn<-fsn(yi = log(incidents_bound_final$RIRR), vi = incidents_bound_final$RIRR_var, type = "Orwin", target = 0.05)
orwin_fsn

### Egger regression
egger_test <- regtest(log(incidents_bound_final$RIRR), vi = incidents_bound_final$RIRR_var, model = 'lm', predictor="sei")
egger_test
########################################
### 7) Forest Plot
########################################
###First round the RES, variance, pooled RES, and pooled variance to 2 decimal places for the foresplot structure
fp.RIRR <- round(RIRR, digits = 2)
fp.lower.CI <- round(lower.CI, digits = 2)
fp.upper.CI <- round(upper.CI, digits = 2)
fp.mva_RIRR <- round(mva_RIRR, digits = 2)
fp.mva_RIRR_CI_low <- round(mva_RIRR_CI_low, digits = 2)
fp.mva_RIRR_CI_upp <- round(mva_RIRR_CI_upp, digits = 2)
fp.z.p.value <- format.pval(p, eps = .001, digits = 2)
fp.mva_RIRR_z_p_value <- format.pval(mva_ts_p, eps = .001, digits = 2)
fp.I2 <- round(I_mva, 2)

###Create the text objects to be displayed (where different from those above). 
fp.txt.RIRR <- format(fp.RIRR, nsmall = 2)
fp.txt.lower.CI <- paste0(format(fp.lower.CI, nsmall = 2), ",")
fp.txt.CI <-  paste0("[", (paste(fp.txt.lower.CI, format(fp.upper.CI, nsmall = 2))), "]")
fp.txt.mva_RIRR_CI_low <- paste0(round(mva_RIRR_CI_low, digits = 2), ",")
fp.txt.mva.CI <-  paste0("[", (paste(fp.txt.mva_RIRR_CI_low , round(mva_RIRR_CI_upp, digits = 2))), "]")

###Amend the year and author observations
fp_authors <- str_remove(incidents_bound_final$authors, ", P. M.,")
fp_author.year <- paste(fp_authors, paste0("(", incidents_bound_final$Year, ")"))

### Create Forest Plot Structure and Text
fp_structure <- 
  structure(list(
    mean  = c(NA, NA, NA, fp.RIRR, NA, fp.mva_RIRR), 
    lower = c(NA, NA, NA, fp.lower.CI, NA, fp.mva_RIRR_CI_low),
    upper = c(NA, NA, NA, fp.upper.CI, NA, fp.mva_RIRR_CI_upp)),
    .Names = c("RES", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame")

fp_tabletext<-cbind(
  c("", "Study", "", fp_author.year , NA, "Summary"),
  c("", "RIRR", "", fp.txt.RIRR, NA, fp.mva_RIRR),
  c("", "95% confidence", "interval", fp.txt.CI, NA, fp.txt.mva.CI),
  c("", "p-Value", "", fp.z.p.value, NA, fp.mva_RIRR_z_p_value))

# trellis.device(device="windows", height = 25, width = 40, color=TRUE)

forestplot(fp_tabletext,
           hrzl_lines = list("4" = gpar(lwd=1, col = "#000044"), 
                             "27" = gpar(lwd=1, col = "#000044"),
                             "28" = gpar(lwd=2, col = "#000044")),
           graph.pos = 4,
           graphwidth = unit(10, "cm"),
           xticks = c(0.25, 0.50, 1.0, 2.0, 4.0, 8.0, 12.0), 
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=.6), ticks = gpar(cex=.6)),
           fp_structure,new_page = TRUE,
           is.summary=c(FALSE, TRUE, TRUE, rep(FALSE,23),TRUE), 
           xlog=TRUE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
           align = c("l", "c", "c", "r"),
           ci.vertices = TRUE)

x <- unit(0.05, 'npc')
y <- unit(.043, 'npc')
# g <- grid.locator('npc') # to select on the plot view. replace x and y with g$x and g$y in grid.text below
grid.text(substitute(paste('Model:', a, '= ',b ), list(a= ~ I^{2}, b=fp.I2)), x, y, just = "left", gp = gpar(fontsize=7, font = 1))
########################################
### 8) Fixed effect model with categorical (Multiplicative Variance Adjustment) and continuous moderators 
########################################
###Create a new data frame for the moderator analysis. 
mods<-incidents_bound_final

#LOG(es)
mods$LRIRR <- log(mods$RIRR)

#Drop levels
mods<- droplevels(mods)

########################################
### 8.1) Methodological
########################################
########################################
### 8.1.1) Year of publication
########################################
###Create a new variable that is numeric and excludes letters from 'Year' 
mods$year_recode <- as.character(mods$Year)
mods$year_recode <- gsub("[^0-9]","", x = mods$year_recode)
mods$year_recode <- as.numeric(mods$year_recode)

### Create an object for the row name, level_ref, and level_alternative names
year_names<- c('year_of_publication', 'methodological', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable', 
               'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
               'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.year<-lm(LRIRR ~ 1 + mods$year_recode,weights = 1/v)

###Extract moderator Beta
mod_year_B_LRIRR<-summary(model.lm.mod.year)$coefficients[2, 1]

###Extract the Beta's standard error
mod_year_B_se<-summary(model.lm.mod.year)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_year_res_df<- model.lm.mod.year$df.residual

###Calculate Qres for the model
mod_year_Y <- as.matrix(LRIRR)
mod_year_W <- diag(1/v)
mod_year_X <- model.matrix(model.lm.mod.year)
mod_year_P <- mod_year_W - mod_year_W %*% mod_year_X %*% solve(t(mod_year_X) %*% mod_year_W %*% mod_year_X) %*% t(mod_year_X) %*% mod_year_W 

mod_year_Qres <- max(0, c(crossprod(mod_year_Y, mod_year_P) %*% mod_year_Y))

###Bind the list of values to the dataframe 'moderator_output'
year_names_list <- c(year_names, mod_year_B_LRIRR, mod_year_B_se, 
                    mod_year_res_df, mod_year_Qres)

names(year_names_list) <-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                          "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v", "alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                          "mod_Qres")
moderator_output <- year_names_list
########################################
### 8.1.2) Retrospective vs Contemporary studygroups
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_studygroups_names<- c('studygroups','methodological', 'retrospective', 'standard') 
### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_studygroups_ret_es <-mods[mods$Studygroups== 'retrospective', "LRIRR"]
sg_studygroups_ret_v <-mods[mods$Studygroups== 'retrospective', "RIRR_var"]
sg_studygroups_ret_k <-length(sg_studygroups_ret_es)
sg_studygroups_ret_PLRIRR <-(sum(sg_studygroups_ret_es/sg_studygroups_ret_v))/(sum(1/sg_studygroups_ret_v))
sg_studygroups_ret_PLRIRR_v <- 1/sum(1/sg_studygroups_ret_v)
sg_studygroups_ret_Q <-sum((sg_studygroups_ret_es-sg_studygroups_ret_PLRIRR)^2/sg_studygroups_ret_v)
sg_studygroups_ret_Q_p <-pchisq(sg_studygroups_ret_Q, df=(length(sg_studygroups_ret_es)-1), lower.tail = FALSE)
sg_studygroups_ret_I2 <-((sg_studygroups_ret_Q-(length(sg_studygroups_ret_es)-1))/sg_studygroups_ret_Q)*100

sg_studygroups_contem_es <-mods[mods$Studygroups== 'standard', "LRIRR"]
sg_studygroups_contem_v <-mods[mods$Studygroups== 'standard', "RIRR_var"]
sg_studygroups_contem_k <-length(sg_studygroups_contem_es)
sg_studygroups_contem_PLRIRR <-(sum(sg_studygroups_contem_es/sg_studygroups_contem_v))/(sum(1/sg_studygroups_contem_v))
sg_studygroups_contem_PLRIRR_v <- 1/sum(1/sg_studygroups_contem_v)
sg_studygroups_contem_Q <-sum((sg_studygroups_contem_es-sg_studygroups_contem_PLRIRR)^2/sg_studygroups_contem_v)
sg_studygroups_contem_Q_p <-pchisq(sg_studygroups_contem_Q, df=(length(sg_studygroups_contem_es)-1), lower.tail = FALSE)
sg_studygroups_contem_I2 <-((sg_studygroups_contem_Q-(length(sg_studygroups_contem_es)-1))/sg_studygroups_contem_Q)*100

###Insert the moderator into the model
model.lm.mod.sour.pub<-lm(LRIRR ~ 1 + mods$Studygroups,weights = 1/v)

###Extract moderator Beta
mod_studygroups_LRIRR<-summary(model.lm.mod.sour.pub)$coefficients[2, 1]

###Extract the Beta's standard error
mod_studygroups_se<-summary(model.lm.mod.sour.pub)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_studygroups_res_df<- model.lm.mod.sour.pub$df.residual

###Calculate Qres for the model
mod_studygroups_Y <- as.matrix(mods$LRIRR)
mod_studygroups_W <- diag(1/mods$RIRR_var)
mod_studygroups_X <- model.matrix(model.lm.mod.sour.pub)
mod_studygroups_P <- mod_studygroups_W - mod_studygroups_W %*% mod_studygroups_X %*% solve(
  t(mod_studygroups_X) %*% mod_studygroups_W %*% mod_studygroups_X) %*% t(
    mod_studygroups_X) %*% mod_studygroups_W 

mod_studygroups_Qres <- max(0, c(crossprod(mod_studygroups_Y, mod_studygroups_P) %*% mod_studygroups_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_studygroups_list<- c(mod_studygroups_names, sg_studygroups_ret_k, sg_studygroups_ret_PLRIRR, sg_studygroups_ret_PLRIRR_v,  
                      sg_studygroups_ret_Q, sg_studygroups_ret_Q_p, sg_studygroups_ret_I2, sg_studygroups_contem_k,
                      sg_studygroups_contem_PLRIRR, sg_studygroups_contem_PLRIRR_v, sg_studygroups_contem_Q, sg_studygroups_contem_Q_p,
                      sg_studygroups_contem_I2, mod_studygroups_LRIRR, mod_studygroups_se, 
                      mod_studygroups_res_df, mod_studygroups_Qres)

names(mod_studygroups_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                            "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                            "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_studygroups_list)
########################################
### 8.1.3) Pre and post period (months)
########################################
###Convert to factor
mods$pt_period_months <- as.factor(mods$pt_period_months)

### Create an object for the row name, level_ref, and level_alternative names
pre_post_period_names<- c('pre_and_post_period_lengths','methodological', '12 months', '24 months') 

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_pre_post_period_12_es <-na.omit(mods[mods$pt_period_months == '12', "LRIRR"])
sg_pre_post_period_12_v <-na.omit(mods[mods$pt_period_months == '12', "RIRR_var"])
sg_pre_post_period_12_k <-length(sg_pre_post_period_12_es)
sg_pre_post_period_12_PLRIRR <-(sum(sg_pre_post_period_12_es/sg_pre_post_period_12_v))/(sum(1/sg_pre_post_period_12_v))
sg_pre_post_period_12_PLRIRR_v <- 1/sum(1/sg_pre_post_period_12_v)
sg_pre_post_period_12_Q <-sum((sg_pre_post_period_12_es-sg_pre_post_period_12_PLRIRR)^2/sg_pre_post_period_12_v)
sg_pre_post_period_12_Q_p <-pchisq(sg_pre_post_period_12_Q, df=(length(sg_pre_post_period_12_es)-1), lower.tail = FALSE)
sg_pre_post_period_12_I2 <-((sg_pre_post_period_12_Q-(length(sg_pre_post_period_12_es)-1))/sg_pre_post_period_12_Q)*100

sg_pre_post_period_24_es <-na.omit(mods[mods$pt_period_months == '24', "LRIRR"])
sg_pre_post_period_24_v <-na.omit(mods[mods$pt_period_months == '24', "RIRR_var"])
sg_pre_post_period_24_k <-length(sg_pre_post_period_24_es)
sg_pre_post_period_24_PLRIRR <-(sum(sg_pre_post_period_24_es/sg_pre_post_period_24_v))/(sum(1/sg_pre_post_period_24_v))
sg_pre_post_period_24_PLRIRR_v <- 1/sum(1/sg_pre_post_period_24_v)
sg_pre_post_period_24_Q <-sum((sg_pre_post_period_24_es-sg_pre_post_period_24_PLRIRR)^2/sg_pre_post_period_24_v)
sg_pre_post_period_24_Q_p <-pchisq(sg_pre_post_period_24_Q, df=(length(sg_pre_post_period_24_es)-1), lower.tail = FALSE)
sg_pre_post_period_24_I2 <-((sg_pre_post_period_24_Q-(length(sg_pre_post_period_24_es)-1))/sg_pre_post_period_24_Q)*100
###Recode 18 months as NA 
mods$pt_period_months <-na_if(mods$pt_period_months, '18')
mods$pt_period_months <-droplevels(mods$pt_period_months)

###Insert the moderator into the model
model.lm.mod.period<-lm(LRIRR ~ 1 + mods$pt_period_months,weights = 1/v)

###Extract moderator Beta
mod_pre_post_period_LRIRR<-summary(model.lm.mod.period)$coefficients[2, 1]

###Extract the Beta's standard error
mod_pre_post_period_se<-summary(model.lm.mod.period)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_pre_post_period_res_df<- model.lm.mod.period$df.residual

###Calculate Qres for the model
mod_pre_post_es <-mods[mods$pt_period_months %in% c("12", "24"), "LRIRR"]
mod_pre_post_v <-mods[mods$pt_period_months %in% c("12", "24"), "RIRR_var"]

mod_pre_post_period_Y <- as.matrix(mod_pre_post_es)
mod_pre_post_period_W <- diag(1/mod_pre_post_v)
mod_pre_post_period_X <- model.matrix(model.lm.mod.period)
mod_pre_post_period_P <- mod_pre_post_period_W - mod_pre_post_period_W %*% mod_pre_post_period_X %*% solve(
  t(mod_pre_post_period_X) %*% mod_pre_post_period_W %*% mod_pre_post_period_X) %*% t(
    mod_pre_post_period_X) %*% mod_pre_post_period_W 

mod_pre_post_period_Qres <- max(0, c(crossprod(mod_pre_post_period_Y, mod_pre_post_period_P) %*% mod_pre_post_period_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_pre_post_period_list<- c(pre_post_period_names,sg_pre_post_period_12_k, sg_pre_post_period_12_PLRIRR,sg_pre_post_period_12_PLRIRR_v, 
                             sg_pre_post_period_12_Q, sg_pre_post_period_12_Q_p, sg_pre_post_period_12_I2, 
                             sg_pre_post_period_24_k, sg_pre_post_period_24_PLRIRR,sg_pre_post_period_24_PLRIRR_v, sg_pre_post_period_24_Q, sg_pre_post_period_24_Q_p,
                             sg_pre_post_period_24_I2, mod_pre_post_period_LRIRR, mod_pre_post_period_se, 
                             mod_pre_post_period_res_df, mod_pre_post_period_Qres)

names(mod_pre_post_period_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                   "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                   "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_pre_post_period_list)
########################################
### 8.1.4) Months measured post court-program
########################################
### Create an object for the row name, level_ref, and level_alternative names
post_prog_names<- c('months_measured_post-program', 'methodological', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable', 
                    'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                    'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.post.prog<-lm(LRIRR ~ 1 + mods$Post_program_PT_months_treat,weights = 1/v)

###Extract moderator Beta
mod_post_prog_B_LRIRR<-summary(model.lm.mod.post.prog)$coefficients[2, 1]

###Extract the Beta's standard error
mod_post_prog_B_se<-summary(model.lm.mod.post.prog)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_post_prog_res_df<- model.lm.mod.post.prog$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$Post_program_PT_months_treat_na <- is.na(mods$Post_program_PT_months_treat)
post_program_PT_months_treat_es <- mods[mods$Post_program_PT_months_treat_na == 'FALSE', "LRIRR"]
post_program_PT_months_treat_v <- mods[mods$Post_program_PT_months_treat_na == 'FALSE', "RIRR_var"]

mod_post_prog_Y <- as.matrix(post_program_PT_months_treat_es)
mod_post_prog_W <- diag(1/post_program_PT_months_treat_v)
mod_post_prog_X <- model.matrix(model.lm.mod.post.prog)
mod_post_prog_P <- mod_post_prog_W - mod_post_prog_W %*% mod_post_prog_X %*% solve(
  t(mod_post_prog_X) %*% mod_post_prog_W %*% mod_post_prog_X) %*% t(
    mod_post_prog_X) %*% mod_post_prog_W 

mod_post_prog_Qres <- max(0, c(crossprod(mod_post_prog_Y, mod_post_prog_P) %*% mod_post_prog_Y))

###Bind the list of values to the dataframe 'moderator_output'
post_prog_list<- c(post_prog_names, mod_post_prog_B_LRIRR, mod_post_prog_B_se, 
                    mod_post_prog_res_df, mod_post_prog_Qres)

names(post_prog_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                          "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                          "mod_Qres")
moderator_output<- bind_rows(moderator_output, post_prog_list)
########################################
### 8.1.5) Published/unpublished
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_sour_pub_names<- c('Source_published','methodological', 'published', 'unpublished') 
### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_sour_pub_pub_es <-mods[mods$source_published== 'published', "LRIRR"]
sg_sour_pub_pub_v <-mods[mods$source_published== 'published', "RIRR_var"]
sg_sour_pub_pub_k <-length(sg_sour_pub_pub_es)
sg_sour_pub_pub_PLRIRR <-(sum(sg_sour_pub_pub_es/sg_sour_pub_pub_v))/(sum(1/sg_sour_pub_pub_v))
sg_sour_pub_pub_PLRIRR_v <- 1/sum(1/sg_sour_pub_pub_v)
sg_sour_pub_pub_Q <-sum((sg_sour_pub_pub_es-sg_sour_pub_pub_PLRIRR)^2/sg_sour_pub_pub_v)
sg_sour_pub_pub_Q_p <-pchisq(sg_sour_pub_pub_Q, df=(length(sg_sour_pub_pub_es)-1), lower.tail = FALSE)
sg_sour_pub_pub_I2 <-((sg_sour_pub_pub_Q-(length(sg_sour_pub_pub_es)-1))/sg_sour_pub_pub_Q)*100

sg_sour_pub_unpub_es <-mods[mods$source_published== 'unpublished', "LRIRR"]
sg_sour_pub_unpub_v <-mods[mods$source_published== 'unpublished', "RIRR_var"]
sg_sour_pub_unpub_k <-length(sg_sour_pub_unpub_es)
sg_sour_pub_unpub_PLRIRR <-(sum(sg_sour_pub_unpub_es/sg_sour_pub_unpub_v))/(sum(1/sg_sour_pub_unpub_v))
sg_sour_pub_unpub_PLRIRR_v <- 1/sum(1/sg_sour_pub_unpub_v)
sg_sour_pub_unpub_Q <-sum((sg_sour_pub_unpub_es-sg_sour_pub_unpub_PLRIRR)^2/sg_sour_pub_unpub_v)
sg_sour_pub_unpub_Q_p <-pchisq(sg_sour_pub_unpub_Q, df=(length(sg_sour_pub_unpub_es)-1), lower.tail = FALSE)
sg_sour_pub_unpub_I2 <-((sg_sour_pub_unpub_Q-(length(sg_sour_pub_unpub_es)-1))/sg_sour_pub_unpub_Q)*100

###Insert the moderator into the model
model.lm.mod.sour.pub<-lm(LRIRR ~ 1 + mods$source_published,weights = 1/v)

###Extract moderator Beta
mod_sour_pub_LRIRR<-summary(model.lm.mod.sour.pub)$coefficients[2, 1]

###Extract the Beta's standard error
mod_sour_pub_se<-summary(model.lm.mod.sour.pub)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_sour_pub_res_df<- model.lm.mod.sour.pub$df.residual

###Calculate Qres for the model
mod_sour_pub_Y <- as.matrix(mods$LRIRR)
mod_sour_pub_W <- diag(1/mods$RIRR_var)
mod_sour_pub_X <- model.matrix(model.lm.mod.sour.pub)
mod_sour_pub_P <- mod_sour_pub_W - mod_sour_pub_W %*% mod_sour_pub_X %*% solve(
  t(mod_sour_pub_X) %*% mod_sour_pub_W %*% mod_sour_pub_X) %*% t(
    mod_sour_pub_X) %*% mod_sour_pub_W 

mod_sour_pub_Qres <- max(0, c(crossprod(mod_sour_pub_Y, mod_sour_pub_P) %*% mod_sour_pub_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_sour_pub_list<- c(mod_sour_pub_names, sg_sour_pub_pub_k, sg_sour_pub_pub_PLRIRR,sg_sour_pub_pub_PLRIRR_v,
                      sg_sour_pub_pub_Q, sg_sour_pub_pub_Q_p, sg_sour_pub_pub_I2, sg_sour_pub_unpub_k,
                      sg_sour_pub_unpub_PLRIRR, sg_sour_pub_unpub_PLRIRR_v, sg_sour_pub_unpub_Q, sg_sour_pub_unpub_Q_p,
                      sg_sour_pub_unpub_I2, mod_sour_pub_LRIRR, mod_sour_pub_se, 
                      mod_sour_pub_res_df, mod_sour_pub_Qres)

names(mod_sour_pub_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                            "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                            "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_sour_pub_list)
########################################
### 8.1.6) Published by NPC Research 
########################################
###Recode variable
mods$organisation_recode <-recode(mods$organisation_or_research_group, "Policy Research Associates" = "None")
# mods$organisation_recode <-droplevels(mods$organisation_recode)

### Create an object for the row name, level_ref, and level_alternative names
mod_org_names<- c('Source_organisation','methodological', 'None', 'NPC Research') 
### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_org_none_es <-mods[mods$organisation_recode== 'None', "LRIRR"]
sg_org_none_v <-mods[mods$organisation_recode== 'None', "RIRR_var"]
sg_org_none_k <-length(sg_org_none_es)
sg_org_none_PLRIRR <-(sum(sg_org_none_es/sg_org_none_v))/(sum(1/sg_org_none_v))
sg_org_none_PLRIRR_v <- 1/sum(1/sg_org_none_v)
sg_org_none_Q <-sum((sg_org_none_es-sg_org_none_PLRIRR)^2/sg_org_none_v)
sg_org_none_Q_p <-pchisq(sg_org_none_Q, df=(length(sg_org_none_es)-1), lower.tail = FALSE)
sg_org_none_I2 <-((sg_org_none_Q-(length(sg_org_none_es)-1))/sg_org_none_Q)*100

sg_org_NPC_es <-mods[mods$organisation_recode== 'NPC research', "LRIRR"]
sg_org_NPC_v <-mods[mods$organisation_recode== 'NPC research', "RIRR_var"]
sg_org_NPC_k <-length(sg_org_NPC_es)
sg_org_NPC_PLRIRR <-(sum(sg_org_NPC_es/sg_org_NPC_v))/(sum(1/sg_org_NPC_v))
sg_org_NPC_PLRIRR_v <- 1/sum(1/sg_org_NPC_v)
sg_org_NPC_Q <-sum((sg_org_NPC_es-sg_org_NPC_PLRIRR)^2/sg_org_NPC_v)
sg_org_NPC_Q_p <-pchisq(sg_org_NPC_Q, df=(length(sg_org_NPC_es)-1), lower.tail = FALSE)
sg_org_NPC_I2 <-((sg_org_NPC_Q-(length(sg_org_NPC_es)-1))/sg_org_NPC_Q)*100

###Insert the moderator into the model
model.lm.mod.org<-lm(LRIRR ~ 1 + mods$organisation_recode,weights = 1/v)

###Extract moderator Beta
mod_org_LRIRR<-summary(model.lm.mod.org)$coefficients[2, 1]

###Extract the Beta's standard error
mod_org_se<-summary(model.lm.mod.org)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_org_res_df<- model.lm.mod.org$df.residual

###Calculate Qres for the model
mod_org_Y <- as.matrix(mods$LRIRR)
mod_org_W <- diag(1/mods$RIRR_var)
mod_org_X <- model.matrix(model.lm.mod.org)
mod_org_P <- mod_org_W - mod_org_W %*% mod_org_X %*% solve(
  t(mod_org_X) %*% mod_org_W %*% mod_org_X) %*% t(
    mod_org_X) %*% mod_org_W 

mod_org_Qres <- max(0, c(crossprod(mod_org_Y, mod_org_P) %*% mod_org_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_org_list<- c(mod_org_names, sg_org_none_k, sg_org_none_PLRIRR, sg_org_none_PLRIRR_v,
                      sg_org_none_Q, sg_org_none_Q_p, sg_org_none_I2, sg_org_NPC_k,
                      sg_org_NPC_PLRIRR,sg_org_NPC_PLRIRR_v, sg_org_NPC_Q, sg_org_NPC_Q_p,
                      sg_org_NPC_I2, mod_org_LRIRR, mod_org_se, 
                      mod_org_res_df, mod_org_Qres)

names(mod_org_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                            "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                            "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_org_list)

########################################
### 8.2) Bias and quality
########################################
########################################
### 8.2.1) Score on the modified Maryland Scientific Methods Scale
########################################
###Convert the MSMS variable to factor
mods$MSMS <- as.factor(mods$MSMS)
### Create an object for the row name, level_ref, and level_alternative names
mod_MSMS_rating_names <-c('MSMS_rating','Quality_and_Bias', '3', '4')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_mod_MSMS_rating_3_es <-mods[mods$MSMS== '3', "LRIRR"]
sg_mod_MSMS_rating_3_v <-mods[mods$MSMS== '3', "RIRR_var"]
sg_mod_MSMS_rating_3_k <-length(sg_mod_MSMS_rating_3_es)
sg_mod_MSMS_rating_3_PLRIRR <-(sum(sg_mod_MSMS_rating_3_es/sg_mod_MSMS_rating_3_v))/(sum(1/sg_mod_MSMS_rating_3_v))
sg_mod_MSMS_rating_3_PLRIRR_v <- 1/sum(1/sg_mod_MSMS_rating_3_v)
sg_mod_MSMS_rating_3_Q <-sum((sg_mod_MSMS_rating_3_es-sg_mod_MSMS_rating_3_PLRIRR)^2/sg_mod_MSMS_rating_3_v)
sg_mod_MSMS_rating_3_Q_p <-pchisq(sg_mod_MSMS_rating_3_Q, df=(length(sg_mod_MSMS_rating_3_es)-1), lower.tail = FALSE)
sg_mod_MSMS_rating_3_I2 <-((sg_mod_MSMS_rating_3_Q-(length(sg_mod_MSMS_rating_3_es)-1))/sg_mod_MSMS_rating_3_Q)*100 

sg_mod_MSMS_rating_4_es <-mods[mods$MSMS== '4', "LRIRR"]
sg_mod_MSMS_rating_4_v <-mods[mods$MSMS== '4', "RIRR_var"]
sg_mod_MSMS_rating_4_k <-length(sg_mod_MSMS_rating_4_es)
sg_mod_MSMS_rating_4_PLRIRR <-(sum(sg_mod_MSMS_rating_4_es/sg_mod_MSMS_rating_4_v))/(sum(1/sg_mod_MSMS_rating_4_v))
sg_mod_MSMS_rating_4_PLRIRR_v <- 1/sum(1/sg_mod_MSMS_rating_4_v)
sg_mod_MSMS_rating_4_Q <-sum((sg_mod_MSMS_rating_4_es-sg_mod_MSMS_rating_4_PLRIRR)^2/sg_mod_MSMS_rating_4_v)
sg_mod_MSMS_rating_4_Q_p <-pchisq(sg_mod_MSMS_rating_4_Q, df=(length(sg_mod_MSMS_rating_4_es)-1), lower.tail = FALSE)
sg_mod_MSMS_rating_4_I2 <-((sg_mod_MSMS_rating_4_Q-(length(sg_mod_MSMS_rating_4_es)-1))/sg_mod_MSMS_rating_4_Q)*100 

###Insert the moderator into the model
model.lm.mod.MSMS<-lm(LRIRR ~ 1 + mods$MSMS, weights = 1/v)

###Extract moderator Beta
mod_MSMS_rating_LRIRR<-summary(model.lm.mod.MSMS)$coefficients[2, 1]

###Extract the Beta's standard error
mod_MSMS_rating_se<-(summary(model.lm.mod.MSMS)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_MSMS_rating_res_df<- model.lm.mod.MSMS$df.residual

###Calculate Qres for the model
mod_MSMS_rating_Y <- as.matrix(mods$LRIRR)
mod_MSMS_rating_W <- diag(1/mods$RIRR_var)
mod_MSMS_rating_X <- model.matrix(model.lm.mod.MSMS)
mod_MSMS_rating_P <- mod_MSMS_rating_W - mod_MSMS_rating_W %*% mod_MSMS_rating_X %*% solve(
  t(mod_MSMS_rating_X) %*% mod_MSMS_rating_W %*% mod_MSMS_rating_X) %*% t(
    mod_MSMS_rating_X) %*% mod_MSMS_rating_W 

mod_MSMS_rating_Qres <- max(0, c(crossprod(mod_MSMS_rating_Y, mod_MSMS_rating_P) %*% mod_MSMS_rating_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_MSMS_rating_list<-c(mod_MSMS_rating_names,sg_mod_MSMS_rating_3_k, sg_mod_MSMS_rating_3_PLRIRR,sg_mod_MSMS_rating_3_PLRIRR_v, sg_mod_MSMS_rating_3_Q, sg_mod_MSMS_rating_3_Q_p, 
                        sg_mod_MSMS_rating_3_I2,sg_mod_MSMS_rating_4_k, sg_mod_MSMS_rating_4_PLRIRR,sg_mod_MSMS_rating_4_PLRIRR_v, sg_mod_MSMS_rating_4_Q, 
                        sg_mod_MSMS_rating_4_Q_p, sg_mod_MSMS_rating_4_I2, mod_MSMS_rating_LRIRR, mod_MSMS_rating_se,
                        mod_MSMS_rating_res_df, mod_MSMS_rating_Qres)

names(mod_MSMS_rating_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                               "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                               "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_MSMS_rating_list)
########################################
### 8.2.2) ROBINS confounding 
########################################
### Note: all outcomes were rated as 'serious' in this domain and so no moderator analysis was performed.
########################################
### 8.2.3) ROBINS selection of participants rating - Moderate and serious 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_ROBIN_partici_names <-c('ROBINS_selection_rating','Quality_and_Bias', 'moderate', 'serious')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_mod_ROB_partici_mod_es <-mods[mods$ROB_participants== 'moderate', "LRIRR"]
sg_mod_ROB_partici_mod_v <-mods[mods$ROB_participants== 'moderate', "RIRR_var"]
sg_mod_ROB_partici_mod_k <-length(sg_mod_ROB_partici_mod_es)
sg_mod_ROB_partici_mod_PLRIRR <-(sum(sg_mod_ROB_partici_mod_es/sg_mod_ROB_partici_mod_v))/(sum(1/sg_mod_ROB_partici_mod_v))
sg_mod_ROB_partici_mod_PLRIRR_v <- 1/sum(1/sg_mod_ROB_partici_mod_v)
sg_mod_ROB_partici_mod_Q <-sum((sg_mod_ROB_partici_mod_es-sg_mod_ROB_partici_mod_PLRIRR)^2/sg_mod_ROB_partici_mod_v)
sg_mod_ROB_partici_mod_Q_p <-pchisq(sg_mod_ROB_partici_mod_Q, df=(length(sg_mod_ROB_partici_mod_es)-1), lower.tail = FALSE)
sg_mod_ROB_partici_mod_I2 <-((sg_mod_ROB_partici_mod_Q-(length(sg_mod_ROB_partici_mod_es)-1))/sg_mod_ROB_partici_mod_Q)*100 

sg_mod_ROB_partici_ser_es <-mods[mods$ROB_participants== 'serious', "LRIRR"]
sg_mod_ROB_partici_ser_v <-mods[mods$ROB_participants== 'serious', "RIRR_var"]
sg_mod_ROB_partici_ser_k <-length(sg_mod_ROB_partici_ser_es)
sg_mod_ROB_partici_ser_PLRIRR <-(sum(sg_mod_ROB_partici_ser_es/sg_mod_ROB_partici_ser_v))/(sum(1/sg_mod_ROB_partici_ser_v))
sg_mod_ROB_partici_ser_PLRIRR_v <- 1/sum(1/sg_mod_ROB_partici_ser_v)
sg_mod_ROB_partici_ser_Q <-sum((sg_mod_ROB_partici_ser_es-sg_mod_ROB_partici_ser_PLRIRR)^2/sg_mod_ROB_partici_ser_v)
sg_mod_ROB_partici_ser_Q_p <-pchisq(sg_mod_ROB_partici_ser_Q, df=(length(sg_mod_ROB_partici_ser_es)-1), lower.tail = FALSE)
sg_mod_ROB_partici_ser_I2 <-((sg_mod_ROB_partici_ser_Q-(length(sg_mod_ROB_partici_ser_es)-1))/sg_mod_ROB_partici_ser_Q)*100 

###Insert the moderator into the model
model.lm.mod.ROB.partic<-lm(LRIRR ~ 1 + mods$ROB_participants, weights = 1/v)

###Extract moderator Beta
mod_ROBIN_partici_LRIRR<-summary(model.lm.mod.ROB.partic)$coefficients[2, 1]

###Extract the Beta's standard error
mod_ROBIN_partici_se<-(summary(model.lm.mod.ROB.partic)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_ROBIN_partici_res_df<- model.lm.mod.ROB.partic$df.residual

###Calculate Qres for the model
mod_ROBIN_partici_Y <- as.matrix(mods$LRIRR)
mod_ROBIN_partici_W <- diag(1/mods$RIRR_var)
mod_ROBIN_partici_X <- model.matrix(model.lm.mod.ROB.partic)
mod_ROBIN_partici_P <- mod_ROBIN_partici_W - mod_ROBIN_partici_W %*% mod_ROBIN_partici_X %*% solve(
  t(mod_ROBIN_partici_X) %*% mod_ROBIN_partici_W %*% mod_ROBIN_partici_X) %*% t(
    mod_ROBIN_partici_X) %*% mod_ROBIN_partici_W 

mod_ROBIN_partici_Qres <- max(0, c(crossprod(mod_ROBIN_partici_Y, mod_ROBIN_partici_P) %*% mod_ROBIN_partici_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_ROBIN_partici_list<-c(mod_ROBIN_partici_names,sg_mod_ROB_partici_mod_k, sg_mod_ROB_partici_mod_PLRIRR,sg_mod_ROB_partici_mod_PLRIRR_v, sg_mod_ROB_partici_mod_Q, sg_mod_ROB_partici_mod_Q_p, 
                          sg_mod_ROB_partici_mod_I2,sg_mod_ROB_partici_ser_k, sg_mod_ROB_partici_ser_PLRIRR,sg_mod_ROB_partici_ser_PLRIRR_v, sg_mod_ROB_partici_ser_Q, 
                          sg_mod_ROB_partici_ser_Q_p, sg_mod_ROB_partici_ser_I2, mod_ROBIN_partici_LRIRR, mod_ROBIN_partici_se,
                          mod_ROBIN_partici_res_df, mod_ROBIN_partici_Qres)

names(mod_ROBIN_partici_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                 "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                 "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_ROBIN_partici_list)

########################################
### 8.2.4) ROBINS classification of treatment status rating - Low and moderate
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_ROBIN_class_names <-c('ROBINS_classification_rating','Quality_and_Bias', 'low', 'moderate')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_mod_ROB_class_low_es <-mods[mods$ROB_classification== 'low', "LRIRR"]
sg_mod_ROB_class_low_v <-mods[mods$ROB_classification== 'low', "RIRR_var"]
sg_mod_ROB_class_low_k <-length(sg_mod_ROB_class_low_es)
sg_mod_ROB_class_low_PLRIRR <-(sum(sg_mod_ROB_class_low_es/sg_mod_ROB_class_low_v))/(sum(1/sg_mod_ROB_class_low_v))
sg_mod_ROB_class_low_PLRIRR_v <- 1/sum(1/sg_mod_ROB_class_low_v)
sg_mod_ROB_class_low_Q <-sum((sg_mod_ROB_class_low_es-sg_mod_ROB_class_low_PLRIRR)^2/sg_mod_ROB_class_low_v)
sg_mod_ROB_class_low_Q_p <-pchisq(sg_mod_ROB_class_low_Q, df=(length(sg_mod_ROB_class_low_es)-1), lower.tail = FALSE)
sg_mod_ROB_class_low_I2 <-((sg_mod_ROB_class_low_Q-(length(sg_mod_ROB_class_low_es)-1))/sg_mod_ROB_class_low_Q)*100 

sg_mod_ROB_class_mod_es <-mods[mods$ROB_classification== 'moderate', "LRIRR"]
sg_mod_ROB_class_mod_v <-mods[mods$ROB_classification== 'moderate', "RIRR_var"]
sg_mod_ROB_class_mod_k <-length(sg_mod_ROB_class_mod_es)
sg_mod_ROB_class_mod_PLRIRR <-(sum(sg_mod_ROB_class_mod_es/sg_mod_ROB_class_mod_v))/(sum(1/sg_mod_ROB_class_mod_v))
sg_mod_ROB_class_mod_PLRIRR_v <- 1/sum(1/sg_mod_ROB_class_mod_v)
sg_mod_ROB_class_mod_Q <-sum((sg_mod_ROB_class_mod_es-sg_mod_ROB_class_mod_PLRIRR)^2/sg_mod_ROB_class_mod_v)
sg_mod_ROB_class_mod_Q_p <-pchisq(sg_mod_ROB_class_mod_Q, df=(length(sg_mod_ROB_class_mod_es)-1), lower.tail = FALSE)
sg_mod_ROB_class_mod_I2 <-((sg_mod_ROB_class_mod_Q-(length(sg_mod_ROB_class_mod_es)-1))/sg_mod_ROB_class_mod_Q)*100 

###Insert the moderator into the model
model.lm.mod.ROB.class<-lm(LRIRR ~ 1 + mods$ROB_classification, weights = 1/v)

###Extract moderator Beta
mod_ROBIN_class_LRIRR<-summary(model.lm.mod.ROB.class)$coefficients[2, 1]

###Extract the Beta's standard error
mod_ROBIN_class_se<-(summary(model.lm.mod.ROB.class)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_ROBIN_class_res_df<- model.lm.mod.ROB.class$df.residual

###Calculate Qres for the model
mod_ROBIN_class_Y <- as.matrix(mods$LRIRR)
mod_ROBIN_class_W <- diag(1/mods$RIRR_var)
mod_ROBIN_class_X <- model.matrix(model.lm.mod.ROB.class)
mod_ROBIN_class_P <- mod_ROBIN_class_W - mod_ROBIN_class_W %*% mod_ROBIN_class_X %*% solve(
  t(mod_ROBIN_class_X) %*% mod_ROBIN_class_W %*% mod_ROBIN_class_X) %*% t(
    mod_ROBIN_class_X) %*% mod_ROBIN_class_W 

mod_ROBIN_class_Qres <- max(0, c(crossprod(mod_ROBIN_class_Y, mod_ROBIN_class_P) %*% mod_ROBIN_class_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_ROBIN_class_list<-c(mod_ROBIN_class_names,sg_mod_ROB_class_mod_k, sg_mod_ROB_class_mod_PLRIRR,sg_mod_ROB_class_low_PLRIRR_v, sg_mod_ROB_class_mod_Q, sg_mod_ROB_class_mod_Q_p, 
                        sg_mod_ROB_class_mod_I2,sg_mod_ROB_class_low_k, sg_mod_ROB_class_low_PLRIRR, sg_mod_ROB_class_mod_PLRIRR_v, sg_mod_ROB_class_low_Q, 
                        sg_mod_ROB_class_low_Q_p, sg_mod_ROB_class_low_I2, mod_ROBIN_class_LRIRR, mod_ROBIN_class_se,
                        mod_ROBIN_class_res_df, mod_ROBIN_class_Qres)

names(mod_ROBIN_class_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                               "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                               "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_ROBIN_class_list)
########################################
### 8.2.5) ROBINS deviation from intended intervention
########################################
### Note: a moderator analysis of the ROBINS deviation ratings was not performed
###       because there was only sufficient information to rate two domains, and they were given two different ratings - 
###       - low and moderate. 
########################################
### 8.2.6.1) ROBINS missing rating - Low and moderate
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_ROBIN_miss_names_moderate <-c('ROBINS_missing_rating','Quality_and_Bias', 'low', 'moderate')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_mod_ROB_miss_low_es <-mods[mods$ROB_missing== 'low', "LRIRR"]
sg_mod_ROB_miss_low_v <-mods[mods$ROB_missing== 'low', "RIRR_var"]
sg_mod_ROB_miss_low_k <-length(sg_mod_ROB_miss_low_es)
sg_mod_ROB_miss_low_PLRIRR <-(sum(sg_mod_ROB_miss_low_es/sg_mod_ROB_miss_low_v))/(sum(1/sg_mod_ROB_miss_low_v))
sg_mod_ROB_miss_low_PLRIRR_v <- 1/sum(1/sg_mod_ROB_miss_low_v)
sg_mod_ROB_miss_low_Q <-sum((sg_mod_ROB_miss_low_es-sg_mod_ROB_miss_low_PLRIRR)^2/sg_mod_ROB_miss_low_v)
sg_mod_ROB_miss_low_Q_p <-pchisq(sg_mod_ROB_miss_low_Q, df=(length(sg_mod_ROB_miss_low_es)-1), lower.tail = FALSE)
sg_mod_ROB_miss_low_I2 <-((sg_mod_ROB_miss_low_Q-(length(sg_mod_ROB_miss_low_es)-1))/sg_mod_ROB_miss_low_Q)*100 

sg_mod_ROB_miss_mod_es <-mods[mods$ROB_missing== 'moderate', "LRIRR"]
sg_mod_ROB_miss_mod_v <-mods[mods$ROB_missing== 'moderate', "RIRR_var"]
sg_mod_ROB_miss_mod_k <-length(sg_mod_ROB_miss_mod_es)
sg_mod_ROB_miss_mod_PLRIRR <-(sum(sg_mod_ROB_miss_mod_es/sg_mod_ROB_miss_mod_v))/(sum(1/sg_mod_ROB_miss_mod_v))
sg_mod_ROB_miss_mod_PLRIRR_v <- 1/sum(1/sg_mod_ROB_miss_mod_v)
sg_mod_ROB_miss_mod_Q <-sum((sg_mod_ROB_miss_mod_es-sg_mod_ROB_miss_mod_PLRIRR)^2/sg_mod_ROB_miss_mod_v)
sg_mod_ROB_miss_mod_Q_p <-pchisq(sg_mod_ROB_miss_mod_Q, df=(length(sg_mod_ROB_miss_mod_es)-1), lower.tail = FALSE)
sg_mod_ROB_miss_mod_I2 <-((sg_mod_ROB_miss_mod_Q-(length(sg_mod_ROB_miss_mod_es)-1))/sg_mod_ROB_miss_mod_Q)*100 

###Insert the moderator into the model
model.lm.mod.ROB.miss<-lm(LRIRR ~ 1 + mods$ROB_missing, weights = 1/v)

###Extract moderator Beta
mod_ROBIN_miss_mod_LRIRR<-summary(model.lm.mod.ROB.miss)$coefficients[2, 1]

###Extract the Beta's standard error
mod_ROBIN_miss_mod_se<-(summary(model.lm.mod.ROB.miss)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_ROBIN_miss_mod_res_df<- model.lm.mod.ROB.miss$df.residual

###Calculate Qres for the model
mod_ROBIN_miss_Y <- as.matrix(mods$LRIRR)
mod_ROBIN_miss_W <- diag(1/mods$RIRR_var)
mod_ROBIN_miss_X <- model.matrix(model.lm.mod.ROB.miss)
mod_ROBIN_miss_P <- mod_ROBIN_miss_W - mod_ROBIN_miss_W %*% mod_ROBIN_miss_X %*% solve(
  t(mod_ROBIN_miss_X) %*% mod_ROBIN_miss_W %*% mod_ROBIN_miss_X) %*% t(
    mod_ROBIN_miss_X) %*% mod_ROBIN_miss_W 

mod_ROBIN_miss_Qres <- max(0, c(crossprod(mod_ROBIN_miss_Y, mod_ROBIN_miss_P) %*% mod_ROBIN_miss_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_ROBIN_miss_list<-c(mod_ROBIN_miss_names_moderate,sg_mod_ROB_miss_mod_k, sg_mod_ROB_miss_mod_PLRIRR,sg_mod_ROB_miss_low_PLRIRR_v, sg_mod_ROB_miss_mod_Q, sg_mod_ROB_miss_mod_Q_p, 
                       sg_mod_ROB_miss_mod_I2,sg_mod_ROB_miss_low_k, sg_mod_ROB_miss_low_PLRIRR,sg_mod_ROB_miss_mod_PLRIRR_v, sg_mod_ROB_miss_low_Q, 
                       sg_mod_ROB_miss_low_Q_p, sg_mod_ROB_miss_low_I2, mod_ROBIN_miss_mod_LRIRR, mod_ROBIN_miss_mod_se,
                       mod_ROBIN_miss_mod_res_df, mod_ROBIN_miss_Qres)

names(mod_ROBIN_miss_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                              "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                              "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_ROBIN_miss_list)
########################################
### 8.2.6.2) ROBINS missing rating - Low and serious
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_ROBIN_miss_names_serious <-c('ROBINS_missing_rating','Quality_and_Bias', 'low', 'serious')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_mod_ROB_miss_low_es <-mods[mods$ROB_missing== 'low', "LRIRR"]
sg_mod_ROB_miss_low_v <-mods[mods$ROB_missing== 'low', "RIRR_var"]
sg_mod_ROB_miss_low_k <-length(sg_mod_ROB_miss_low_es)
sg_mod_ROB_miss_low_PLRIRR <-(sum(sg_mod_ROB_miss_low_es/sg_mod_ROB_miss_low_v))/(sum(1/sg_mod_ROB_miss_low_v))
sg_mod_ROB_miss_low_PLRIRR_v <- 1/sum(1/sg_mod_ROB_miss_low_v)
sg_mod_ROB_miss_low_Q <-sum((sg_mod_ROB_miss_low_es-sg_mod_ROB_miss_low_PLRIRR)^2/sg_mod_ROB_miss_low_v)
sg_mod_ROB_miss_low_Q_p <-pchisq(sg_mod_ROB_miss_low_Q, df=(length(sg_mod_ROB_miss_low_es)-1), lower.tail = FALSE)
sg_mod_ROB_miss_low_I2 <-((sg_mod_ROB_miss_low_Q-(length(sg_mod_ROB_miss_low_es)-1))/sg_mod_ROB_miss_low_Q)*100 

sg_mod_ROB_miss_ser_es <-mods[mods$ROB_missing== 'serious', "LRIRR"]
sg_mod_ROB_miss_ser_v <-mods[mods$ROB_missing== 'serious', "RIRR_var"]
sg_mod_ROB_miss_ser_k <-length(sg_mod_ROB_miss_ser_es)
sg_mod_ROB_miss_ser_PLRIRR <-(sum(sg_mod_ROB_miss_ser_es/sg_mod_ROB_miss_ser_v))/(sum(1/sg_mod_ROB_miss_ser_v))
sg_mod_ROB_miss_ser_PLRIRR_v <- 1/sum(1/sg_mod_ROB_miss_ser_v)
sg_mod_ROB_miss_ser_Q <-sum((sg_mod_ROB_miss_ser_es-sg_mod_ROB_miss_ser_PLRIRR)^2/sg_mod_ROB_miss_ser_v)
sg_mod_ROB_miss_ser_Q_p <-pchisq(sg_mod_ROB_miss_ser_Q, df=(length(sg_mod_ROB_miss_ser_es)-1), lower.tail = FALSE)
sg_mod_ROB_miss_ser_I2 <-((sg_mod_ROB_miss_ser_Q-(length(sg_mod_ROB_miss_ser_es)-1))/sg_mod_ROB_miss_ser_Q)*100 

###Insert the moderator into the model
model.lm.mod.ROB.miss<-lm(LRIRR ~ 1 + mods$ROB_missing, weights = 1/v)

###Extract moderator Beta
mod_ROBIN_miss_ser_LRIRR<-summary(model.lm.mod.ROB.miss)$coefficients[3, 1]

###Extract the Beta's standard error
mod_ROBIN_miss_ser_se<-(summary(model.lm.mod.ROB.miss)$coefficients[3, 2])

###Create an object for the residual degrees of freedom. 
mod_ROBIN_miss_ser_df<- model.lm.mod.ROB.miss$df.residual

###Calculate Qres for the model
mod_ROBIN_miss_Y <- as.matrix(mods$LRIRR)
mod_ROBIN_miss_W <- diag(1/mods$RIRR_var)
mod_ROBIN_miss_X <- model.matrix(model.lm.mod.ROB.miss)
mod_ROBIN_miss_P <- mod_ROBIN_miss_W - mod_ROBIN_miss_W %*% mod_ROBIN_miss_X %*% solve(
  t(mod_ROBIN_miss_X) %*% mod_ROBIN_miss_W %*% mod_ROBIN_miss_X) %*% t(
    mod_ROBIN_miss_X) %*% mod_ROBIN_miss_W 

mod_ROBIN_miss_Qres <- max(0, c(crossprod(mod_ROBIN_miss_Y, mod_ROBIN_miss_P) %*% mod_ROBIN_miss_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_ROBIN_miss_list<-c(mod_ROBIN_miss_names_serious,sg_mod_ROB_miss_ser_k, sg_mod_ROB_miss_ser_PLRIRR,sg_mod_ROB_miss_low_PLRIRR_v, sg_mod_ROB_miss_ser_Q, sg_mod_ROB_miss_ser_Q_p, 
                       sg_mod_ROB_miss_ser_I2,sg_mod_ROB_miss_low_k, sg_mod_ROB_miss_low_PLRIRR,sg_mod_ROB_miss_ser_PLRIRR_v, sg_mod_ROB_miss_low_Q, 
                       sg_mod_ROB_miss_low_Q_p, sg_mod_ROB_miss_low_I2, mod_ROBIN_miss_ser_LRIRR, mod_ROBIN_miss_ser_se,
                       mod_ROBIN_miss_ser_df, mod_ROBIN_miss_Qres)

names(mod_ROBIN_miss_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                              "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                              "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_ROBIN_miss_list)
########################################
### 8.2.7) ROBINS measurement rating - Low and moderate
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_ROBIN_meas_names <-c('ROBINS_measurement_rating','Quality_and_Bias', 'low', 'moderate')

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_mod_ROB_meas_low_es <-na.omit(mods[mods$ROB_measurement== 'low', "LRIRR"])
sg_mod_ROB_meas_low_v <-na.omit(mods[mods$ROB_measurement== 'low', "RIRR_var"])
sg_mod_ROB_meas_low_k <-length(sg_mod_ROB_meas_low_es)
sg_mod_ROB_meas_low_PLRIRR <-(sum(sg_mod_ROB_meas_low_es/sg_mod_ROB_meas_low_v))/(sum(1/sg_mod_ROB_meas_low_v))
sg_mod_ROB_meas_low_PLRIRR_v <- 1/sum(1/sg_mod_ROB_meas_low_v)
sg_mod_ROB_meas_low_Q <-sum((sg_mod_ROB_meas_low_es-sg_mod_ROB_meas_low_PLRIRR)^2/sg_mod_ROB_meas_low_v)
sg_mod_ROB_meas_low_Q_p <-pchisq(sg_mod_ROB_meas_low_Q, df=(length(sg_mod_ROB_meas_low_es)-1), lower.tail = FALSE)
sg_mod_ROB_meas_low_I2 <-((sg_mod_ROB_meas_low_Q-(length(sg_mod_ROB_meas_low_es)-1))/sg_mod_ROB_meas_low_Q)*100 

sg_mod_ROB_meas_mod_es <-na.omit(mods[mods$ROB_measurement== 'moderate', "LRIRR"])
sg_mod_ROB_meas_mod_v <-na.omit(mods[mods$ROB_measurement== 'moderate', "RIRR_var"])
sg_mod_ROB_meas_mod_k <-length(sg_mod_ROB_meas_mod_es)
sg_mod_ROB_meas_mod_PLRIRR <-(sum(sg_mod_ROB_meas_mod_es/sg_mod_ROB_meas_mod_v))/(sum(1/sg_mod_ROB_meas_mod_v))
sg_mod_ROB_meas_mod_PLRIRR_v <- 1/sum(1/sg_mod_ROB_meas_mod_v)
sg_mod_ROB_meas_mod_Q <-sum((sg_mod_ROB_meas_mod_es-sg_mod_ROB_meas_mod_PLRIRR)^2/sg_mod_ROB_meas_mod_v)
sg_mod_ROB_meas_mod_Q_p <-pchisq(sg_mod_ROB_meas_mod_Q, df=(length(sg_mod_ROB_meas_mod_es)-1), lower.tail = FALSE)
sg_mod_ROB_meas_mod_I2 <-((sg_mod_ROB_meas_mod_Q-(length(sg_mod_ROB_meas_mod_es)-1))/sg_mod_ROB_meas_mod_Q)*100 

###Recode NI months as NA 
mods$ROB_measurement <-na_if(mods$ROB_measurement, 'NI')
# mods$ROB_measurement <-droplevels(mods$ROB_measurement)

###Insert the moderator into the model
model.lm.mod.ROB.meas<-lm(LRIRR ~ 1 + mods$ROB_measurement, weights = 1/v)

###Extract moderator Beta
mod_ROBIN_meas_LRIRR<-summary(model.lm.mod.ROB.meas)$coefficients[2, 1]

###Extract the Beta's standard error
mod_ROBIN_meas_se<-summary(model.lm.mod.ROB.meas)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_ROBIN_meas_res_df<- model.lm.mod.ROB.meas$df.residual

###Calculate Qres for the model
mod_ROBIN_meas_es <-mods[mods$ROB_measurement %in% c("low", "moderate"), "LRIRR"]
mod_ROBIN_meas_v <-mods[mods$ROB_measurement %in% c("low", "moderate"), "RIRR_var"]

mod_ROBIN_meas_Y <- as.matrix(mod_ROBIN_meas_es)
mod_ROBIN_meas_W <- diag(1/mod_ROBIN_meas_v)
mod_ROBIN_meas_X <- model.matrix(model.lm.mod.ROB.meas)
mod_ROBIN_meas_P <- mod_ROBIN_meas_W - mod_ROBIN_meas_W %*% mod_ROBIN_meas_X %*% solve(
  t(mod_ROBIN_meas_X) %*% mod_ROBIN_meas_W %*% mod_ROBIN_meas_X) %*% t(
    mod_ROBIN_meas_X) %*% mod_ROBIN_meas_W 

mod_ROBIN_meas_Qres <- max(0, c(crossprod(mod_ROBIN_meas_Y, mod_ROBIN_meas_P) %*% mod_ROBIN_meas_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_ROBIN_meas_list<- c(mod_ROBIN_meas_names,sg_mod_ROB_meas_low_k, sg_mod_ROB_meas_low_PLRIRR, sg_mod_ROB_meas_low_PLRIRR_v,
                        sg_mod_ROB_meas_low_Q, sg_mod_ROB_meas_low_Q_p, sg_mod_ROB_meas_low_I2, 
                        sg_mod_ROB_meas_mod_k, sg_mod_ROB_meas_mod_PLRIRR, sg_mod_ROB_meas_mod_PLRIRR_v, sg_mod_ROB_meas_mod_Q, sg_mod_ROB_meas_mod_Q_p,
                        sg_mod_ROB_meas_mod_I2, mod_ROBIN_meas_LRIRR, mod_ROBIN_meas_se, 
                        mod_ROBIN_meas_res_df, mod_ROBIN_meas_Qres)

names(mod_ROBIN_meas_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                              "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                              "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_ROBIN_meas_list)
########################################
### 8.2.8) ROBINS selection of reported result 
########################################
### Note: all outcomes were rated as 'moderate' in this domain and so no moderator analysis was performed
########################################
### 8.3) Sample characteristics  
########################################
########################################
### 8.3.1) Lifetime arrests - Treatment 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_lt_arr_tr_names<- c('lifetime_arrests_treat', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                        'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                        'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.lt.arr.tr<-lm(LRIRR ~ 1 + mods$LT_arrests_treat, weights = 1/v)

###Extract moderator Beta
mod_lt_arr_tr_B_LRIRR<-summary(model.lm.mod.lt.arr.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_lt_arr_tr_B_se<-summary(model.lm.mod.lt.arr.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_lt_arr_tr_res_df<- model.lm.mod.lt.arr.tr$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$LT_arrests_treat_na<- is.na(mods$LT_arrests_treat)
mod_lt_arr_tr_es <- mods[mods$LT_arrests_treat_na == 'FALSE', "LRIRR"]
mod_lt_arr_tr_v <- mods[mods$LT_arrests_treat_na == 'FALSE', "RIRR_var"]

mod_lt_arr_tr_Y <- as.matrix(mod_lt_arr_tr_es)
mod_lt_arr_tr_W <- diag(1/mod_lt_arr_tr_v)
mod_lt_arr_tr_X <- model.matrix(model.lm.mod.lt.arr.tr)
mod_lt_arr_tr_P <- mod_lt_arr_tr_W - mod_lt_arr_tr_W %*% mod_lt_arr_tr_X %*% solve(
  t(mod_lt_arr_tr_X) %*% mod_lt_arr_tr_W %*% mod_lt_arr_tr_X) %*% t(
    mod_lt_arr_tr_X) %*% mod_lt_arr_tr_W 

mod_lt_arr_tr_Qres <- max(0, c(crossprod(mod_lt_arr_tr_Y, mod_lt_arr_tr_P) %*% mod_lt_arr_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_lt_arr_tr_names_list<- c(mod_lt_arr_tr_names, mod_lt_arr_tr_B_LRIRR, mod_lt_arr_tr_B_se, 
                             mod_lt_arr_tr_res_df, mod_lt_arr_tr_Qres)

names(mod_lt_arr_tr_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                   "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                   "mod_Qres")
moderator_output<- bind_rows(moderator_output, mod_lt_arr_tr_names_list)
########################################
### 8.3.2) Lifetime arrests - Control
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_lt_arr_con_names<- c('lifetime_arrests_cont', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                         'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                         'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.lt.arr.con<-lm(LRIRR ~ 1 + mods$LT_arrests_cont, weights = 1/v)

###Extract moderator Beta
mod_lt_arr_con_B_LRIRR<-summary(model.lm.mod.lt.arr.con)$coefficients[2, 1]

###Extract the Beta's standard error
mod_lt_arr_con_B_se<-summary(model.lm.mod.lt.arr.con)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_lt_arr_con_res_df<- model.lm.mod.lt.arr.con$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$LT_arrests_cont_na<- is.na(mods$LT_arrests_cont)
mod_lt_arr_con_es <- mods[mods$LT_arrests_cont_na == 'FALSE', "LRIRR"]
mod_lt_arr_con_v <- mods[mods$LT_arrests_cont_na == 'FALSE', "RIRR_var"]

mod_lt_arr_con_Y <- as.matrix(mod_lt_arr_con_es)
mod_lt_arr_con_W <- diag(1/mod_lt_arr_con_v)
mod_lt_arr_con_X <- model.matrix(model.lm.mod.lt.arr.con)
mod_lt_arr_con_P <- mod_lt_arr_con_W - mod_lt_arr_con_W %*% mod_lt_arr_con_X %*% solve(
  t(mod_lt_arr_con_X) %*% mod_lt_arr_con_W %*% mod_lt_arr_con_X) %*% t(
    mod_lt_arr_con_X) %*% mod_lt_arr_con_W 

mod_lt_arr_con_Qres <- max(0, c(crossprod(mod_lt_arr_con_Y, mod_lt_arr_con_P) %*% mod_lt_arr_con_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_lt_arr_con_names_list<- c(mod_lt_arr_con_names, mod_lt_arr_con_B_LRIRR, mod_lt_arr_con_B_se, 
                              mod_lt_arr_con_res_df, mod_lt_arr_con_Qres)

names(mod_lt_arr_con_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                    "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                    "mod_Qres")
moderator_output<- bind_rows(moderator_output, mod_lt_arr_con_names_list)
########################################
### 8.3.3) Mean age - Treatment 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_m_age_tr_names<- c('mean_age_treat', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                       'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                       'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.m.age.tr<-lm(LRIRR ~ 1 + mods$mean_age_treat, weights = 1/v)

###Extract moderator Beta
mod_m_age_tr_B_LRIRR<-summary(model.lm.mod.m.age.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_m_age_tr_B_se<-summary(model.lm.mod.m.age.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_m_age_tr_res_df<- model.lm.mod.m.age.tr$df.residual

###Calculate Qres for the model
mod_m_age_tr_Y <- as.matrix(LRIRR)
mod_m_age_tr_W <- diag(1/v)
mod_m_age_tr_X <- model.matrix(model.lm.mod.m.age.tr)
mod_m_age_tr_P <- mod_m_age_tr_W - mod_m_age_tr_W %*% mod_m_age_tr_X %*% solve(
  t(mod_m_age_tr_X) %*% mod_m_age_tr_W %*% mod_m_age_tr_X) %*% t(
    mod_m_age_tr_X) %*% mod_m_age_tr_W 

mod_m_age_tr_Qres <- max(0, c(crossprod(mod_m_age_tr_Y, mod_m_age_tr_P) %*% mod_m_age_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_m_age_tr_names_list<- c(mod_m_age_tr_names, mod_m_age_tr_B_LRIRR, mod_m_age_tr_B_se, 
                            mod_m_age_tr_res_df, mod_m_age_tr_Qres)

names(mod_m_age_tr_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                  "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                  "mod_Qres")
moderator_output<- bind_rows(moderator_output, mod_m_age_tr_names_list)
########################################
### 8.3.4) Mean age - Control
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_m_age_con_names<- c('mean_age_cont', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                        'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                        'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.m.age.con<-lm(LRIRR ~ 1 + mods$mean_age_cont, weights = 1/v)

###Extract moderator Beta
mod_m_age_con_B_LRIRR<-summary(model.lm.mod.m.age.con)$coefficients[2, 1]

###Extract the Beta's standard error
mod_m_age_con_B_se<-summary(model.lm.mod.m.age.con)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_m_age_con_res_df<- model.lm.mod.m.age.con$df.residual

###Calculate Qres for the model
mod_m_age_con_Y <- as.matrix(LRIRR)
mod_m_age_con_W <- diag(1/v)
mod_m_age_con_X <- model.matrix(model.lm.mod.m.age.con)
mod_m_age_con_P <- mod_m_age_con_W - mod_m_age_con_W %*% mod_m_age_con_X %*% solve(
  t(mod_m_age_con_X) %*% mod_m_age_con_W %*% mod_m_age_con_X) %*% t(
    mod_m_age_con_X) %*% mod_m_age_con_W 

mod_m_age_con_Qres <- max(0, c(crossprod(mod_m_age_con_Y, mod_m_age_con_P) %*% mod_m_age_con_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_m_age_con_names_list<- c(mod_m_age_con_names, mod_m_age_con_B_LRIRR, mod_m_age_con_B_se, 
                             mod_m_age_con_res_df, mod_m_age_con_Qres)

names(mod_m_age_con_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                   "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                   "mod_Qres")
moderator_output<- bind_rows(moderator_output, mod_m_age_con_names_list)
########################################
### 8.3.5) Percent male - Treatment 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_percent_male_tr_names<- c('percent_male_treat', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                              'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                              'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.percent.male.tr<-lm(LRIRR ~ 1 + mods$male_._treat, weights = 1/v)

###Extract moderator Beta
mod_percent_male_tr_B_LRIRR<-summary(model.lm.mod.percent.male.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_percent_male_tr_B_se<-summary(model.lm.mod.percent.male.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_percent_male_tr_res_df<- model.lm.mod.percent.male.tr$df.residual

###Calculate Qres for the model
mod_percent_male_tr_Y <- as.matrix(LRIRR)
mod_percent_male_tr_W <- diag(1/v)
mod_percent_male_tr_X <- model.matrix(model.lm.mod.percent.male.tr)
mod_percent_male_tr_P <- mod_percent_male_tr_W - mod_percent_male_tr_W %*% mod_percent_male_tr_X %*% solve(
  t(mod_percent_male_tr_X) %*% mod_percent_male_tr_W %*% mod_percent_male_tr_X) %*% t(
    mod_percent_male_tr_X) %*% mod_percent_male_tr_W 

mod_percent_male_tr_Qres <- max(0, c(crossprod(mod_percent_male_tr_Y, mod_percent_male_tr_P) %*% mod_percent_male_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_percent_male_tr_names_list<- c(mod_percent_male_tr_names, mod_percent_male_tr_B_LRIRR, mod_percent_male_tr_B_se, 
                                   mod_percent_male_tr_res_df, mod_percent_male_tr_Qres)

names(mod_percent_male_tr_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                         "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                         "mod_Qres")
moderator_output<- bind_rows(moderator_output, mod_percent_male_tr_names_list)
########################################
### 8.3.6) Percent male - Control 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_percent_male_con_names<- c('percent_male_cont', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                               'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                               'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.percent.male.con<-lm(LRIRR ~ 1 + mods$male_._cont, weights = 1/v)

###Extract moderator Beta
mod_percent_male_con_B_LRIRR<-summary(model.lm.mod.percent.male.con)$coefficients[2, 1]

###Extract the Beta's standard error
mod_percent_male_con_B_se<-summary(model.lm.mod.percent.male.con)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_percent_male_con_res_df<- model.lm.mod.percent.male.con$df.residual

###Calculate Qres for the model
mod_percent_male_con_Y <- as.matrix(LRIRR)
mod_percent_male_con_W <- diag(1/v)
mod_percent_male_con_X <- model.matrix(model.lm.mod.percent.male.con)
mod_percent_male_con_P <- mod_percent_male_con_W - mod_percent_male_con_W %*% mod_percent_male_con_X %*% solve(
  t(mod_percent_male_con_X) %*% mod_percent_male_con_W %*% mod_percent_male_con_X) %*% t(
    mod_percent_male_con_X) %*% mod_percent_male_con_W 

mod_percent_male_con_Qres <- max(0, c(crossprod(mod_percent_male_con_Y, mod_percent_male_con_P) %*% mod_percent_male_con_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_percent_male_con_names_list<- c(mod_percent_male_con_names, mod_percent_male_con_B_LRIRR, mod_percent_male_con_B_se, 
                                    mod_percent_male_con_res_df, mod_percent_male_con_Qres)

names(mod_percent_male_con_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                          "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                          "mod_Qres")
moderator_output<- bind_rows(moderator_output, mod_percent_male_con_names_list)
########################################
### 8.3.7) Percent white - Treatment 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_percent_white_tr_names<- c('percent_white_treat', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                               'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                               'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.percent.white.tr<-lm(LRIRR ~ 1 + mods$white_._treat, weights = 1/v)

###Extract moderator Beta
mod_percent_white_tr_B_LRIRR<-summary(model.lm.mod.percent.white.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_percent_white_tr_B_se<-summary(model.lm.mod.percent.white.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_percent_white_tr_res_df<- model.lm.mod.percent.white.tr$df.residual

###Calculate Qres for the model
mod_percent_white_tr_Y <- as.matrix(LRIRR)
mod_percent_white_tr_W <- diag(1/v)
mod_percent_white_tr_X <- model.matrix(model.lm.mod.percent.white.tr)
mod_percent_white_tr_P <- mod_percent_white_tr_W - mod_percent_white_tr_W %*% mod_percent_white_tr_X %*% solve(
  t(mod_percent_white_tr_X) %*% mod_percent_white_tr_W %*% mod_percent_white_tr_X) %*% t(
    mod_percent_white_tr_X) %*% mod_percent_white_tr_W 

mod_percent_white_tr_Qres <- max(0, c(crossprod(mod_percent_white_tr_Y, mod_percent_white_tr_P) %*% mod_percent_white_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_percent_white_tr_names_list<- c(mod_percent_white_tr_names, mod_percent_white_tr_B_LRIRR, mod_percent_white_tr_B_se, 
                                    mod_percent_white_tr_res_df, mod_percent_white_tr_Qres)

names(mod_percent_white_tr_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                          "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                          "mod_Qres")
moderator_output<- bind_rows(moderator_output, mod_percent_white_tr_names_list)
########################################
### 8.3.8) Percent white - Control 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_percent_white_con_names<- c('percent_white_cont', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                                'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                                'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.percent.white.con<-lm(LRIRR ~ 1 + mods$white_._cont, weights = 1/v)

###Extract moderator Beta
mod_percent_white_con_B_LRIRR<-summary(model.lm.mod.percent.white.con)$coefficients[2, 1]

###Extract the Beta's standard error
mod_percent_white_con_B_se<-summary(model.lm.mod.percent.white.con)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_percent_white_con_res_df<- model.lm.mod.percent.white.con$df.residual

###Calculate Qres for the model
mod_percent_white_con_Y <- as.matrix(LRIRR)
mod_percent_white_con_W <- diag(1/v)
mod_percent_white_con_X <- model.matrix(model.lm.mod.percent.white.con)
mod_percent_white_con_P <- mod_percent_white_con_W - mod_percent_white_con_W %*% mod_percent_white_con_X %*% solve(
  t(mod_percent_white_con_X) %*% mod_percent_white_con_W %*% mod_percent_white_con_X) %*% t(
    mod_percent_white_con_X) %*% mod_percent_white_con_W 

mod_percent_white_con_Qres <- max(0, c(crossprod(mod_percent_white_con_Y, mod_percent_white_con_P) %*% mod_percent_white_con_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_percent_white_con_names_list<- c(mod_percent_white_con_names, mod_percent_white_con_B_LRIRR, mod_percent_white_con_B_se, 
                                     mod_percent_white_con_res_df, mod_percent_white_con_Qres)

names(mod_percent_white_con_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                           "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                           "mod_Qres")
moderator_output<- bind_rows(moderator_output, mod_percent_white_con_names_list)
########################################
### 8.3.9) Percent program graduates - Treatment 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_percent_pro_grad_tr_names<- c('percent_program_graduate_treat', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                                  'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                                  'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.percent.prog.grad.tr<-lm(LRIRR ~ 1 + mods$pt_._program_graduates, weights = 1/v)

###Extract moderator Beta
mod_percent_pro_grad_tr_B_LRIRR<-summary(model.lm.mod.percent.prog.grad.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_percent_pro_grad_tr_B_se<-summary(model.lm.mod.percent.prog.grad.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_percent_pro_grad_tr_res_df<- model.lm.mod.percent.prog.grad.tr$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$pt_prog_grad_tr_na<- is.na(mods$pt_._program_graduates)
mod_prog_grad_tr_es <- mods[mods$pt_prog_grad_tr_na == 'FALSE', "LRIRR"]
mod_prog_grad_tr_v <- mods[mods$pt_prog_grad_tr_na == 'FALSE', "RIRR_var"]

mod_percent_pro_grad_tr_Y <- as.matrix(mod_prog_grad_tr_es)
mod_percent_pro_grad_tr_W <- diag(1/mod_prog_grad_tr_v)
mod_percent_pro_grad_tr_X <- model.matrix(model.lm.mod.percent.prog.grad.tr)
mod_percent_pro_grad_tr_P <- mod_percent_pro_grad_tr_W - mod_percent_pro_grad_tr_W %*% mod_percent_pro_grad_tr_X %*% solve(
  t(mod_percent_pro_grad_tr_X) %*% mod_percent_pro_grad_tr_W %*% mod_percent_pro_grad_tr_X) %*% t(
    mod_percent_pro_grad_tr_X) %*% mod_percent_pro_grad_tr_W 

mod_percent_pro_grad_tr_Qres <- max(0, c(crossprod(mod_percent_pro_grad_tr_Y, mod_percent_pro_grad_tr_P) %*% mod_percent_pro_grad_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_percent_pro_grad_tr_list<- c(mod_percent_pro_grad_tr_names, mod_percent_pro_grad_tr_B_LRIRR, mod_percent_pro_grad_tr_B_se, 
                                 mod_percent_pro_grad_tr_res_df, mod_percent_pro_grad_tr_Qres)

names(mod_percent_pro_grad_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                       "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                       "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_percent_pro_grad_tr_list)
########################################
### 8.3.10) Percent program terminates - Treatment 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_percent_pro_term_tr_names<- c('percent_program_terminates_treat', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                                  'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                                  'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.percent.prog.term.tr<-lm(LRIRR ~ 1 + mods$pt_._program_terminates, weights = 1/v)

###Extract moderator Beta
mod_percent_pro_term_tr_B_LRIRR<-summary(model.lm.mod.percent.prog.term.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_percent_pro_term_tr_B_se<-summary(model.lm.mod.percent.prog.term.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_percent_pro_term_tr_res_df<- model.lm.mod.percent.prog.term.tr$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$pt_prog_term_tr_na<- is.na(mods$pt_._program_terminates)
mod_prog_term_tr_es <- mods[mods$pt_prog_term_tr_na == 'FALSE', "LRIRR"]
mod_prog_term_tr_v <- mods[mods$pt_prog_term_tr_na == 'FALSE', "RIRR_var"]

mod_percent_pro_term_tr_Y <- as.matrix(mod_prog_term_tr_es)
mod_percent_pro_term_tr_W <- diag(1/mod_prog_term_tr_v)
mod_percent_pro_term_tr_X <- model.matrix(model.lm.mod.percent.prog.term.tr)
mod_percent_pro_term_tr_P <- mod_percent_pro_term_tr_W - mod_percent_pro_term_tr_W %*% mod_percent_pro_term_tr_X %*% solve(
  t(mod_percent_pro_term_tr_X) %*% mod_percent_pro_term_tr_W %*% mod_percent_pro_term_tr_X) %*% t(
    mod_percent_pro_term_tr_X) %*% mod_percent_pro_term_tr_W 

mod_percent_pro_term_tr_Qres <- max(0, c(crossprod(mod_percent_pro_term_tr_Y, mod_percent_pro_term_tr_P) %*% mod_percent_pro_term_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_percent_pro_term_tr_list<- c(mod_percent_pro_term_tr_names, mod_percent_pro_term_tr_B_LRIRR, mod_percent_pro_term_tr_B_se, 
                                 mod_percent_pro_term_tr_res_df, mod_percent_pro_term_tr_Qres)

names(mod_percent_pro_term_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                       "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                       "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_percent_pro_term_tr_list)
########################################
### 8.3.11) Percent treatment completed - Treatment 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_treat_comp_tr_names<- c('treatment_completion_treat', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                            'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                            'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.treat.comp.tr<-lm(LRIRR ~ 1 + mods$pt_treatment_completion_treat_., weights = 1/v)

###Extract moderator Beta
mod_treat_comp_tr_B_LRIRR<-summary(model.lm.mod.treat.comp.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_treat_comp_tr_B_se<-summary(model.lm.mod.treat.comp.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_treat_comp_tr_res_df<- model.lm.mod.treat.comp.tr$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$pt_treatment_completion_treat_na <- is.na(mods$pt_treatment_completion_treat_.)
mod_treat_comp_tr_es <- mods[mods$pt_treatment_completion_treat_na == 'FALSE', "LRIRR"]
mod_treat_comp_tr_v <- mods[mods$pt_treatment_completion_treat_na == 'FALSE', "RIRR_var"]

mod_treat_comp_tr_Y <- as.matrix(mod_treat_comp_tr_es)
mod_treat_comp_tr_W <- diag(1/mod_treat_comp_tr_v)
mod_treat_comp_tr_X <- model.matrix(model.lm.mod.treat.comp.tr)
mod_treat_comp_tr_P <- mod_treat_comp_tr_W - mod_treat_comp_tr_W %*% mod_treat_comp_tr_X %*% solve(
  t(mod_treat_comp_tr_X) %*% mod_treat_comp_tr_W %*% mod_treat_comp_tr_X) %*% t(
    mod_treat_comp_tr_X) %*% mod_treat_comp_tr_W 

mod_treat_comp_tr_Qres <- max(0, c(crossprod(mod_treat_comp_tr_Y, mod_treat_comp_tr_P) %*% mod_treat_comp_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_treat_comp_tr_list<- c(mod_treat_comp_tr_names, mod_treat_comp_tr_B_LRIRR, mod_treat_comp_tr_B_se, 
                           mod_treat_comp_tr_res_df, mod_treat_comp_tr_Qres)

names(mod_treat_comp_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                 "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                 "mod_Qres")
moderator_output<- bind_rows(moderator_output, mod_treat_comp_tr_list)
########################################
### 8.3.12) Percent treatment completed - Control 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_treat_comp_con_names<- c('treatment_completion_cont', 'sample_characteristics', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                             'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                             'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.treat.comp.tr<-lm(LRIRR ~ 1 + mods$pt_treatment_completion_cont_., weights = 1/v)

###Extract moderator Beta
mod_treat_comp_con_B_LRIRR<-summary(model.lm.mod.treat.comp.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_treat_comp_con_B_se<-summary(model.lm.mod.treat.comp.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_treat_comp_con_res_df<- model.lm.mod.treat.comp.tr$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$pt_treatment_completion_cont_na <- is.na(mods$pt_treatment_completion_cont_.)
mod_treat_comp_con_es <- mods[mods$pt_treatment_completion_cont_na == 'FALSE', "LRIRR"]
mod_treat_comp_con_v <- mods[mods$pt_treatment_completion_cont_na == 'FALSE', "RIRR_var"]

mod_treat_comp_con_Y <- as.matrix(mod_treat_comp_con_es)
mod_treat_comp_con_W <- diag(1/mod_treat_comp_con_v)
mod_treat_comp_con_X <- model.matrix(model.lm.mod.treat.comp.tr)
mod_treat_comp_con_P <- mod_treat_comp_con_W - mod_treat_comp_con_W %*% mod_treat_comp_con_X %*% solve(
  t(mod_treat_comp_con_X) %*% mod_treat_comp_con_W %*% mod_treat_comp_con_X) %*% t(
    mod_treat_comp_con_X) %*% mod_treat_comp_con_W 

mod_treat_comp_con_Qres <- max(0, c(crossprod(mod_treat_comp_con_Y, mod_treat_comp_con_P) %*% mod_treat_comp_con_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_treat_comp_con_list<- c(mod_treat_comp_con_names, mod_treat_comp_con_B_LRIRR, mod_treat_comp_con_B_se, 
                            mod_treat_comp_con_res_df, mod_treat_comp_con_Qres)

names(mod_treat_comp_con_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                  "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                  "mod_Qres")
moderator_output<- bind_rows(moderator_output, mod_treat_comp_con_list)

########################################
### 8.4) Programmatic 
########################################
########################################
### 8.4.1.1) Court type - Drug court
########################################
### If factor, re-level the setting variable for the model
if (class(mods$setting) == 'factor') {relevel(mods$setting, ref = "MHC")}

### Create an object for the row name, level_ref, and level_alternative names
setting_mod_DC_names<- c('court_type', 'programmatic', 'MHC', 'DC')

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_setting_MHC_es <-mods[mods$setting == 'MHC', "LRIRR"]
sg_setting_MHC_v <-mods[mods$setting == 'MHC', "RIRR_var"]
sg_setting_MHC_k <-length(sg_setting_MHC_es)
sg_setting_MHC_PLRIRR <-(sum(sg_setting_MHC_es/sg_setting_MHC_v))/(sum(1/sg_setting_MHC_v))
sg_setting_MHC_PLRIRR_v <- 1/sum(1/sg_setting_MHC_v)
sg_setting_MHC_Q <-sum((sg_setting_MHC_es-sg_setting_MHC_PLRIRR)^2/sg_setting_MHC_v)
sg_setting_MHC_Q_p <-pchisq(sg_setting_MHC_Q, df=(length(sg_setting_MHC_es)-1), lower.tail = FALSE)
sg_setting_MHC_I2 <-((sg_setting_MHC_Q-(length(sg_setting_MHC_es)-1))/sg_setting_MHC_Q)*100

sg_setting_DC_es <-mods[mods$setting == 'DC', "LRIRR"]
sg_setting_DC_v <-mods[mods$setting == 'DC', "RIRR_var"]
sg_setting_DC_k <-length(sg_setting_DC_es)
sg_setting_DC_PLRIRR <-(sum(sg_setting_DC_es/sg_setting_DC_v))/(sum(1/sg_setting_DC_v))
sg_setting_DC_PLRIRR_v <- 1/sum(1/sg_setting_DC_v)
sg_setting_DC_Q <-sum((sg_setting_DC_es-sg_setting_DC_PLRIRR)^2/sg_setting_DC_v)
sg_setting_DC_Q_p <-pchisq(sg_setting_DC_Q, df=(length(sg_setting_DC_es)-1), lower.tail = FALSE)
sg_setting_DC_I2 <-((sg_setting_DC_Q-(length(sg_setting_DC_es)-1))/sg_setting_DC_Q)*100

###Insert the moderator into the model
model.lm.mod.setting<-lm(LRIRR ~ 1 + mods$setting,weights = 1/v)

###Extract moderator Beta
mod_setting_DC_LRIRR<-summary(model.lm.mod.setting)$coefficients[2, 1]

###Extract the Beta's standard error of the model
mod_setting_DC_se<-(summary(model.lm.mod.setting)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_setting_DC_res_df<- model.lm.mod.setting$df.residual

###Calculate Qres for the model
mod_setting_Y <- as.matrix(mods$LRIRR)
mod_setting_W <- diag(1/v)
mod_setting_X <- model.matrix(model.lm.mod.setting)
mod_setting_P <- mod_setting_W - mod_setting_W %*% mod_setting_X %*% solve(
  t(mod_setting_X) %*% mod_setting_W %*% mod_setting_X) %*% t(
    mod_setting_X) %*% mod_setting_W 

mod_setting_Qres <- max(0, c(crossprod(mod_setting_Y, mod_setting_P) %*% mod_setting_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_setting_DC_list<-c(setting_mod_DC_names,sg_setting_MHC_k, sg_setting_MHC_PLRIRR,sg_setting_MHC_PLRIRR_v, sg_setting_MHC_Q, sg_setting_MHC_Q_p, 
                       sg_setting_MHC_I2,sg_setting_DC_k, sg_setting_DC_PLRIRR,sg_setting_DC_PLRIRR_v, sg_setting_DC_Q, sg_setting_DC_Q_p, 
                       sg_setting_DC_I2, mod_setting_DC_LRIRR, mod_setting_DC_se, mod_setting_DC_res_df, 
                       mod_setting_Qres)

names(mod_setting_DC_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                              "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                              "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_setting_DC_list)
########################################
### 8.4.1.2) Court type - Driving-under-the-influence court
########################################
### Create an object for the row name, level_ref, and level_alternative names
setting_mod_DUIC_names<- c('court_type', 'programmatic', 'MHC', 'DUIC')

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroupsg_setting_MHC_es <-mods[mods$setting == 'MHC', "LRIRR"]
sg_setting_MHC_es <-mods[mods$setting == 'MHC', "LRIRR"]
sg_setting_MHC_v <-mods[mods$setting == 'MHC', "RIRR_var"]
sg_setting_MHC_k <-length(sg_setting_MHC_es)
sg_setting_MHC_PLRIRR <-(sum(sg_setting_MHC_es/sg_setting_MHC_v))/(sum(1/sg_setting_MHC_v))
sg_setting_MHC_PLRIRR_v <- 1/sum(1/sg_setting_MHC_v)
sg_setting_MHC_Q <-sum((sg_setting_MHC_es-sg_setting_MHC_PLRIRR)^2/sg_setting_MHC_v)
sg_setting_MHC_Q_p <-pchisq(sg_setting_MHC_Q, df=(length(sg_setting_MHC_es)-1), lower.tail = FALSE)
sg_setting_MHC_I2 <-((sg_setting_MHC_Q-(length(sg_setting_MHC_es)-1))/sg_setting_MHC_Q)*100

sg_setting_DUIC_es <-mods[mods$setting == 'DUIC', "LRIRR"]
sg_setting_DUIC_v <-mods[mods$setting == 'DUIC', "RIRR_var"]
sg_setting_DUIC_k <-length(sg_setting_DUIC_es)
sg_setting_DUIC_PLRIRR <-(sum(sg_setting_DUIC_es/sg_setting_DUIC_v))/(sum(1/sg_setting_DUIC_v))
sg_setting_DUIC_PLRIRR_v <- 1/sum(1/sg_setting_DUIC_v)
sg_setting_DUIC_Q <-sum((sg_setting_DUIC_es-sg_setting_DUIC_PLRIRR)^2/sg_setting_DUIC_v)
sg_setting_DUIC_Q_p <-pchisq(sg_setting_DUIC_Q, df=(length(sg_setting_DUIC_es)-1), lower.tail = FALSE)
sg_setting_DUIC_I2 <-((sg_setting_DUIC_Q-(length(sg_setting_DUIC_es)-1))/sg_setting_DUIC_Q)*100

###Insert the moderator into the model
model.lm.mod.setting<-lm(LRIRR ~ 1 + mods$setting,weights = 1/v)

###Extract moderator Beta
mod_setting_DUIC_LRIRR<-summary(model.lm.mod.setting)$coefficients[3, 1]

###Extract the Beta's standard error
mod_setting_DUIC_se<-(summary(model.lm.mod.setting)$coefficients[3, 2])

###Create an object for the residual degrees of freedom. 
mod_setting_DUIC_res_df<- model.lm.mod.setting$df.residual

###Calculate Qres for the model
mod_setting_Y <- as.matrix(mods$LRIRR)
mod_setting_W <- diag(1/v)
mod_setting_X <- model.matrix(model.lm.mod.setting)
mod_setting_P <- mod_setting_W - mod_setting_W %*% mod_setting_X %*% solve(
  t(mod_setting_X) %*% mod_setting_W %*% mod_setting_X) %*% t(
    mod_setting_X) %*% mod_setting_W 

mod_setting_Qres <- max(0, c(crossprod(mod_setting_Y, mod_setting_P) %*% mod_setting_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_setting_DUIC_list<-c(setting_mod_DUIC_names,sg_setting_MHC_k, sg_setting_MHC_PLRIRR,sg_setting_MHC_PLRIRR_v,  sg_setting_MHC_Q, sg_setting_MHC_Q_p, 
                         sg_setting_MHC_I2,sg_setting_DUIC_k , sg_setting_DUIC_PLRIRR,sg_setting_DUIC_PLRIRR_v, sg_setting_DUIC_Q, sg_setting_DUIC_Q_p, 
                         sg_setting_DUIC_I2, mod_setting_DUIC_LRIRR, mod_setting_DUIC_se, mod_setting_DUIC_res_df, 
                         mod_setting_Qres)

names(mod_setting_DUIC_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_setting_DUIC_list)
########################################
### 8.4.1.3) Court type - Family drug court
########################################
### Create an object for the row name, level_ref, and level_alternative names
setting_mod_FDC_names<- c('court_type','programmatic', 'MHC', 'FDC')

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_setting_MHC_es <-mods[mods$setting == 'MHC', "LRIRR"]
sg_setting_MHC_v <-mods[mods$setting == 'MHC', "RIRR_var"]
sg_setting_MHC_k <-length(sg_setting_MHC_es)
sg_setting_MHC_PLRIRR <-(sum(sg_setting_MHC_es/sg_setting_MHC_v))/(sum(1/sg_setting_MHC_v))
sg_setting_MHC_PLRIRR_v <- 1/sum(1/sg_setting_MHC_v)
sg_setting_MHC_Q <-sum((sg_setting_MHC_es-sg_setting_MHC_PLRIRR)^2/sg_setting_MHC_v)
sg_setting_MHC_Q_p <-pchisq(sg_setting_MHC_Q, df=(length(sg_setting_MHC_es)-1), lower.tail = FALSE)
sg_setting_MHC_I2 <-((sg_setting_MHC_Q-(length(sg_setting_MHC_es)-1))/sg_setting_MHC_Q)*100

sg_setting_FDC_es <-mods[mods$setting == 'Family DC', "LRIRR"]
sg_setting_FDC_v <-mods[mods$setting == 'Family DC', "RIRR_var"]
sg_setting_FDC_k <-length(sg_setting_FDC_es)
sg_setting_FDC_PLRIRR <-(sum(sg_setting_FDC_es/sg_setting_FDC_v))/(sum(1/sg_setting_FDC_v))
sg_setting_FDC_PLRIRR_v <- 1/sum(1/sg_setting_FDC_v)
sg_setting_FDC_Q <-sum((sg_setting_FDC_es-sg_setting_FDC_PLRIRR)^2/sg_setting_FDC_v)
sg_setting_FDC_Q_p <-pchisq(sg_setting_FDC_Q, df=(length(sg_setting_FDC_es)-1), lower.tail = FALSE)
sg_setting_FDC_I2 <-((sg_setting_FDC_Q-(length(sg_setting_FDC_es)-1))/sg_setting_FDC_Q)*100

###Insert the moderator into the model
model.lm.mod.setting<-lm(LRIRR ~ 1 + mods$setting,weights = 1/v)

###Extract moderator Beta
mod_setting_FDC_LRIRR<-summary(model.lm.mod.setting)$coefficients[4, 1]

###Extract the Beta's standard error
mod_setting_FDC_se<-(summary(model.lm.mod.setting)$coefficients[4, 2])

###Create an object for the residual degrees of freedom. 
mod_setting_FDC_res_df<- model.lm.mod.setting$df.residual

###Calculate Qres for the model
mod_setting_Y <- as.matrix(mods$LRIRR)
mod_setting_W <- diag(1/v)
mod_setting_X <- model.matrix(model.lm.mod.setting)
mod_setting_P <- mod_setting_W - mod_setting_W %*% mod_setting_X %*% solve(
  t(mod_setting_X) %*% mod_setting_W %*% mod_setting_X) %*% t(
    mod_setting_X) %*% mod_setting_W 

mod_setting_Qres <- max(0, c(crossprod(mod_setting_Y, mod_setting_P) %*% mod_setting_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_setting_FDC_list<-c(setting_mod_FDC_names, sg_setting_MHC_k,  sg_setting_MHC_PLRIRR,sg_setting_MHC_PLRIRR_v, sg_setting_MHC_Q, sg_setting_MHC_Q_p, 
                        sg_setting_MHC_I2,sg_setting_FDC_k, sg_setting_FDC_PLRIRR,sg_setting_FDC_PLRIRR_v, sg_setting_FDC_Q, sg_setting_FDC_Q_p, 
                        sg_setting_FDC_I2, mod_setting_FDC_LRIRR, mod_setting_FDC_se, mod_setting_FDC_res_df, 
                        mod_setting_Qres)

names(mod_setting_FDC_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                               "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                               "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_setting_FDC_list)
########################################
### 8.4.2) Minimum program duration (months) 
########################################
### Create an object for the row name, level_ref, and level_alternative names
min_durat_names<- c('minimum_program_duration', 'programmatic', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                    'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                    'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.min.durat<-lm(LRIRR ~ 1 + mods$Program_duration_minimum_months, weights = 1/v)

###Extract moderator Beta
mod_min_durat_B_LRIRR<-summary(model.lm.mod.min.durat)$coefficients[2, 1]

###Extract the Beta's standard error
mod_min_durat_B_se<-summary(model.lm.mod.min.durat)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_min_durat_res_df<- model.lm.mod.min.durat$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$Program_duration_minimum_months_na <- is.na(mods$Program_duration_minimum_months)
mod_min_durat_es <- mods[mods$Program_duration_minimum_months_na == 'FALSE', "LRIRR"]
mod_min_durat_v <- mods[mods$Program_duration_minimum_months_na == 'FALSE', "RIRR_var"]

mod_min_durat_Y <- as.matrix(mod_min_durat_es)
mod_min_durat_W <- diag(1/mod_min_durat_v)
mod_min_durat_X <- model.matrix(model.lm.mod.min.durat)
mod_min_durat_P <- mod_min_durat_W - mod_min_durat_W %*% mod_min_durat_X %*% solve(
  t(mod_min_durat_X) %*% mod_min_durat_W %*% mod_min_durat_X) %*% t(
    mod_min_durat_X) %*% mod_min_durat_W 

mod_min_durat_Qres <- max(0, c(crossprod(mod_min_durat_Y, mod_min_durat_P) %*% mod_min_durat_Y))

###Bind the list of values to the dataframe 'moderator_output'
min_durat_list<- c(min_durat_names, mod_min_durat_B_LRIRR, mod_min_durat_B_se, 
                    mod_min_durat_res_df, mod_min_durat_Qres)

names(min_durat_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                          "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                          "mod_Qres")
moderator_output<- bind_rows(moderator_output, min_durat_list)
########################################
### 8.4.3) Average program duration (months) 
########################################
### Create an object for the row name, level_ref, and level_alternative names
avg_durat_names<- c('average_program_duration', 'programmatic', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                    'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                    'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.avg.durat<-lm(LRIRR ~ 1 + mods$Program_duration_average_months,weights = 1/v)

###Extract moderator Beta
mod_avg_durat_B_LRIRR<-summary(model.lm.mod.avg.durat)$coefficients[2, 1]

###Extract the Beta's standard error
mod_avg_durat_B_se<-summary(model.lm.mod.avg.durat)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_avg_durat_res_df<- model.lm.mod.avg.durat$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$Program_duration_average_months_na <- is.na(mods$Program_duration_average_months)
mod_avg_durat_es <- mods[mods$Program_duration_average_months_na == 'FALSE', "LRIRR"]
mod_avg_durat_v <- mods[mods$Program_duration_average_months_na == 'FALSE', "RIRR_var"]

mod_avg_durat_Y <- as.matrix(mod_avg_durat_es)
mod_avg_durat_W <- diag(1/mod_avg_durat_v)
mod_avg_durat_X <- model.matrix(model.lm.mod.avg.durat)
mod_avg_durat_P <- mod_avg_durat_W - mod_avg_durat_W %*% mod_avg_durat_X %*% solve(
  t(mod_avg_durat_X) %*% mod_avg_durat_W %*% mod_avg_durat_X) %*% t(
    mod_avg_durat_X) %*% mod_avg_durat_W 

mod_avg_durat_Qres <- max(0, c(crossprod(mod_avg_durat_Y, mod_avg_durat_P) %*% mod_avg_durat_Y))

###Bind the list of values to the dataframe 'moderator_output'
avg_durat_list<- c(avg_durat_names, mod_avg_durat_B_LRIRR, mod_avg_durat_B_se, 
                    mod_avg_durat_res_df, mod_avg_durat_Qres)

names(avg_durat_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                          "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                          "mod_Qres")

moderator_output<- bind_rows(moderator_output, avg_durat_list)

########################################
### 8.4.4) Accepted violent offenders  
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_acc_vio_names<- c('accepted_violent_offenders','programmatic', 'No', 'Yes')

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_acc_vio_No_es <-na.omit(mods[mods$accepted_violent== 'No', "LRIRR"])
sg_acc_vio_No_v <-na.omit(mods[mods$accepted_violent == 'No', "RIRR_var"])
sg_acc_vio_No_k <-length(sg_acc_vio_No_es)
sg_acc_vio_No_PLRIRR <-(sum(sg_acc_vio_No_es/sg_acc_vio_No_v))/(sum(1/sg_acc_vio_No_v))
sg_acc_vio_No_PLRIRR_v <- 1/sum(1/sg_acc_vio_No_v)
sg_acc_vio_No_Q <-sum((sg_acc_vio_No_es-sg_acc_vio_No_PLRIRR)^2/sg_acc_vio_No_v)
sg_acc_vio_No_Q_p <-pchisq(sg_acc_vio_No_Q, df=(length(sg_acc_vio_No_es)-1), lower.tail = FALSE)
sg_acc_vio_No_I2 <-((sg_acc_vio_No_Q-(length(sg_acc_vio_No_es)-1))/sg_acc_vio_No_Q)*100 

sg_acc_vio_Yes_es <-na.omit(mods[mods$accepted_violent == 'Yes', "LRIRR"])
sg_acc_vio_Yes_v <-na.omit(mods[mods$accepted_violent == 'Yes', "RIRR_var"])
sg_acc_vio_Yes_k <-length(sg_acc_vio_Yes_es)
sg_acc_vio_Yes_PLRIRR <-(sum(sg_acc_vio_Yes_es/sg_acc_vio_Yes_v))/(sum(1/sg_acc_vio_Yes_v))
sg_acc_vio_Yes_PLRIRR_v <- 1/sum(1/sg_acc_vio_Yes_v)
sg_acc_vio_Yes_Q <-sum((sg_acc_vio_Yes_es-sg_acc_vio_Yes_PLRIRR)^2/sg_acc_vio_Yes_v)
sg_acc_vio_Yes_Q_p <-pchisq(sg_acc_vio_Yes_Q, df=(length(sg_acc_vio_Yes_es)-1), lower.tail = FALSE)
sg_acc_vio_Yes_I2 <-((sg_acc_vio_Yes_Q-(length(sg_acc_vio_Yes_es)-1))/sg_acc_vio_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.acc.vio<-lm(LRIRR ~ 1 + mods$accepted_violent,weights = 1/v)

###Extract moderator Beta
mod_acc_vio_LRIRR<-summary(model.lm.mod.acc.vio)$coefficients[2, 1]

###Extract the Beta's standard error
mod_acc_vio_se<-(summary(model.lm.mod.acc.vio)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_acc_vio_res_df<- model.lm.mod.acc.vio$df.residual

###Calculate Qres for the model
mod_acc_vio_es <-mods[mods$accepted_violent %in% c("No", "Yes"), "LRIRR"]
mod_acc_vio_v <-mods[mods$accepted_violent %in% c("No", "Yes"), "RIRR_var"]

mod_acc_vio_Y <- as.matrix(mod_acc_vio_es)
mod_acc_vio_W <- diag(1/mod_acc_vio_v)
mod_acc_vio_X <- model.matrix(model.lm.mod.acc.vio)
mod_acc_vio_P <- mod_acc_vio_W - mod_acc_vio_W %*% mod_acc_vio_X %*% solve(
  t(mod_acc_vio_X) %*% mod_acc_vio_W %*% mod_acc_vio_X) %*% t(
    mod_acc_vio_X) %*% mod_acc_vio_W 

mod_acc_vio_Qres <- max(0, c(crossprod(mod_acc_vio_Y, mod_acc_vio_P) %*% mod_acc_vio_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_acc_vio_list<-c(mod_acc_vio_names,sg_acc_vio_No_k, sg_acc_vio_No_PLRIRR,sg_acc_vio_No_PLRIRR_v, sg_acc_vio_No_Q, sg_acc_vio_No_Q_p, 
                    sg_acc_vio_No_I2,sg_acc_vio_Yes_k, sg_acc_vio_Yes_PLRIRR,sg_acc_vio_Yes_PLRIRR_v, sg_acc_vio_Yes_Q, sg_acc_vio_Yes_Q_p, 
                    sg_acc_vio_Yes_I2, mod_acc_vio_LRIRR, mod_acc_vio_se, mod_acc_vio_res_df, 
                    mod_acc_vio_Qres)

names(mod_acc_vio_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                           "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                           "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_acc_vio_list)
########################################
### 8.4.5) Reported accepting violent offenders  
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_acc_vio_names<- c('reported_accepting_violent_offenders', 'programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_acc_vio_No_es <-mods[mods$reported_accepting_violent== 'No', "LRIRR"]
sg_rep_acc_vio_No_v <-mods[mods$reported_accepting_violent == 'No', "RIRR_var"]
sg_rep_acc_vio_No_k <-length(sg_rep_acc_vio_No_es)
sg_rep_acc_vio_No_PLRIRR <-(sum(sg_rep_acc_vio_No_es/sg_rep_acc_vio_No_v))/(sum(1/sg_rep_acc_vio_No_v))
sg_rep_acc_vio_No_PLRIRR_v <- 1/sum(1/sg_rep_acc_vio_No_v)
sg_rep_acc_vio_No_Q <-sum((sg_rep_acc_vio_No_es-sg_rep_acc_vio_No_PLRIRR)^2/sg_rep_acc_vio_No_v)
sg_rep_acc_vio_No_Q_p <-pchisq(sg_rep_acc_vio_No_Q, df=(length(sg_rep_acc_vio_No_es)-1), lower.tail = FALSE)
sg_rep_acc_vio_No_I2 <-((sg_rep_acc_vio_No_Q-(length(sg_rep_acc_vio_No_es)-1))/sg_rep_acc_vio_No_Q)*100 

sg_rep_acc_vio_Yes_es<- mods[mods$reported_accepting_violent== 'Yes', "LRIRR"]
sg_rep_acc_vio_Yes_v<- mods[mods$reported_accepting_violent == 'Yes', "RIRR_var"]
sg_rep_acc_vio_Yes_k <-length(sg_rep_acc_vio_Yes_es)
sg_rep_acc_vio_Yes_PLRIRR <-(sum(sg_rep_acc_vio_Yes_es/sg_rep_acc_vio_Yes_v))/(sum(1/sg_rep_acc_vio_Yes_v))
sg_rep_acc_vio_Yes_PLRIRR_v <- 1/sum(1/sg_rep_acc_vio_Yes_v)
sg_rep_acc_vio_Yes_Q <-sum((sg_rep_acc_vio_Yes_es-sg_rep_acc_vio_Yes_PLRIRR)^2/sg_rep_acc_vio_Yes_v)
sg_rep_acc_vio_Yes_Q_p<- pchisq(sg_rep_acc_vio_Yes_Q, df=(length(sg_rep_acc_vio_Yes_es)-1), lower.tail = FALSE)
sg_rep_acc_vio_Yes_I2 <-((sg_rep_acc_vio_Yes_Q-(length(sg_rep_acc_vio_Yes_es)-1))/sg_rep_acc_vio_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.acc.vio<-lm(LRIRR ~ 1 + mods$reported_accepting_violent,weights = 1/v)

###Extract moderator Beta
mod_rep_acc_vio_LRIRR<-summary(model.lm.mod.rep.acc.vio)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_acc_vio_se<-(summary(model.lm.mod.rep.acc.vio)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_acc_vio_res_df<- model.lm.mod.rep.acc.vio$df.residual

###Calculate Qres for the model
mod_rep_acc_vio_Y <- as.matrix(mods$LRIRR)
mod_rep_acc_vio_W <- diag(1/mods$RIRR_var)
mod_rep_acc_vio_X <- model.matrix(model.lm.mod.rep.acc.vio)
mod_rep_acc_vio_P <- mod_rep_acc_vio_W - mod_rep_acc_vio_W %*% mod_rep_acc_vio_X %*% solve(
  t(mod_rep_acc_vio_X) %*% mod_rep_acc_vio_W %*% mod_rep_acc_vio_X) %*% t(
    mod_rep_acc_vio_X) %*% mod_rep_acc_vio_W 

mod_rep_acc_vio_Qres <- max(0, c(crossprod(mod_rep_acc_vio_Y, mod_rep_acc_vio_P) %*% mod_rep_acc_vio_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_acc_vio_list<-c(mod_rep_acc_vio_names,sg_rep_acc_vio_No_k, sg_rep_acc_vio_No_PLRIRR,sg_rep_acc_vio_No_PLRIRR_v, sg_rep_acc_vio_No_Q, sg_rep_acc_vio_No_Q_p, 
                        sg_rep_acc_vio_No_I2,sg_rep_acc_vio_Yes_k, sg_rep_acc_vio_Yes_PLRIRR,sg_rep_acc_vio_Yes_PLRIRR_v, sg_rep_acc_vio_Yes_Q, sg_rep_acc_vio_Yes_Q_p, 
                        sg_rep_acc_vio_Yes_I2, mod_rep_acc_vio_LRIRR, mod_rep_acc_vio_se, mod_rep_acc_vio_res_df, mod_rep_acc_vio_Qres)

names(mod_rep_acc_vio_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                               "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                               "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_acc_vio_list)
########################################
### 8.4.6) Accepted felons
########################################
###Note there are insufficient cell numbers to run this moderator analysis.
########################################
### 8.4.7) Reported accepting Felons  
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_acc_fel_names<- c('reported_accepting_felons', 'programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_acc_fel_No_es <-mods[mods$reported_accepting_felons== 'No', "LRIRR"]
sg_rep_acc_fel_No_v <-mods[mods$reported_accepting_felons == 'No', "RIRR_var"]
sg_rep_acc_fel_No_k <-length(sg_rep_acc_fel_No_es)
sg_rep_acc_fel_No_PLRIRR <-(sum(sg_rep_acc_fel_No_es/sg_rep_acc_fel_No_v))/(sum(1/sg_rep_acc_fel_No_v))
sg_rep_acc_fel_No_PLRIRR_v <- 1/sum(1/sg_rep_acc_fel_No_v)
sg_rep_acc_fel_No_Q <-sum((sg_rep_acc_fel_No_es-sg_rep_acc_fel_No_PLRIRR)^2/sg_rep_acc_fel_No_v)
sg_rep_acc_fel_No_Q_p <-pchisq(sg_rep_acc_fel_No_Q, df=(length(sg_rep_acc_fel_No_es)-1), lower.tail = FALSE)
sg_rep_acc_fel_No_I2 <-((sg_rep_acc_fel_No_Q-(length(sg_rep_acc_fel_No_es)-1))/sg_rep_acc_fel_No_Q)*100 

sg_rep_acc_fel_Yes_es<- mods[mods$reported_accepting_felons== 'Yes', "LRIRR"]
sg_rep_acc_fel_Yes_v<- mods[mods$reported_accepting_felons == 'Yes', "RIRR_var"]
sg_rep_acc_fel_Yes_k <-length(sg_rep_acc_fel_Yes_es)
sg_rep_acc_fel_Yes_PLRIRR <-(sum(sg_rep_acc_fel_Yes_es/sg_rep_acc_fel_Yes_v))/(sum(1/sg_rep_acc_fel_Yes_v))
sg_rep_acc_fel_Yes_PLRIRR_v <- 1/sum(1/sg_rep_acc_fel_Yes_v)
sg_rep_acc_fel_Yes_Q <-sum((sg_rep_acc_fel_Yes_es-sg_rep_acc_fel_Yes_PLRIRR)^2/sg_rep_acc_fel_Yes_v)
sg_rep_acc_fel_Yes_Q_p<- pchisq(sg_rep_acc_fel_Yes_Q, df=(length(sg_rep_acc_fel_Yes_es)-1), lower.tail = FALSE)
sg_rep_acc_fel_Yes_I2 <-((sg_rep_acc_fel_Yes_Q-(length(sg_rep_acc_fel_Yes_es)-1))/sg_rep_acc_fel_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.acc.fel<-lm(LRIRR ~ 1 + mods$reported_accepting_felons,weights = 1/v)

###Extract moderator Beta
mod_rep_acc_fel_LRIRR<-summary(model.lm.mod.rep.acc.fel)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_acc_fel_se<-(summary(model.lm.mod.rep.acc.fel)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_acc_fel_res_df<- model.lm.mod.rep.acc.fel$df.residual

###Calculate Qres for the model
mod_rep_acc_fel_Y <- as.matrix(mods$LRIRR)
mod_rep_acc_fel_W <- diag(1/mods$RIRR_var)
mod_rep_acc_fel_X <- model.matrix(model.lm.mod.rep.acc.fel)
mod_rep_acc_fel_P <- mod_rep_acc_fel_W - mod_rep_acc_fel_W %*% mod_rep_acc_fel_X %*% solve(
  t(mod_rep_acc_fel_X) %*% mod_rep_acc_fel_W %*% mod_rep_acc_fel_X) %*% t(
    mod_rep_acc_fel_X) %*% mod_rep_acc_fel_W 

mod_rep_acc_fel_Qres <- max(0, c(crossprod(mod_rep_acc_fel_Y, mod_rep_acc_fel_P) %*% mod_rep_acc_fel_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_acc_fel_list<-c(mod_rep_acc_fel_names,sg_rep_acc_fel_No_k, sg_rep_acc_fel_No_PLRIRR,sg_rep_acc_fel_No_PLRIRR_v, sg_rep_acc_fel_No_Q, sg_rep_acc_fel_No_Q_p, 
                        sg_rep_acc_fel_No_I2,sg_rep_acc_fel_Yes_k, sg_rep_acc_fel_Yes_PLRIRR,sg_rep_acc_fel_Yes_PLRIRR_v, sg_rep_acc_fel_Yes_Q, 
                        sg_rep_acc_fel_Yes_Q_p, sg_rep_acc_fel_Yes_I2, mod_rep_acc_fel_LRIRR, mod_rep_acc_fel_se,
                        mod_rep_acc_fel_res_df, mod_rep_acc_fel_Qres)

names(mod_rep_acc_fel_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                               "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                               "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_acc_fel_list)
########################################
### 8.4.8) Reported that charges were expunged 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_charge_ex_names<- c('reported_expunging_charges','programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_charge_ex_No_es <-mods[mods$reported_charges_expunged== 'No', "LRIRR"]
sg_rep_charge_ex_No_v <-mods[mods$reported_charges_expunged== 'No', "RIRR_var"]
sg_rep_charge_ex_No_k <-length(sg_rep_charge_ex_No_es)
sg_rep_charge_ex_No_PLRIRR <-(sum(sg_rep_charge_ex_No_es/sg_rep_charge_ex_No_v))/(sum(1/sg_rep_charge_ex_No_v))
sg_rep_charge_ex_No_PLRIRR_v <- 1/sum(1/sg_rep_charge_ex_No_v)
sg_rep_charge_ex_No_Q <-sum((sg_rep_charge_ex_No_es-sg_rep_charge_ex_No_PLRIRR)^2/sg_rep_charge_ex_No_v)
sg_rep_charge_ex_No_Q_p <-pchisq(sg_rep_charge_ex_No_Q, df=(length(sg_rep_charge_ex_No_es)-1), lower.tail = FALSE)
sg_rep_charge_ex_No_I2 <-((sg_rep_charge_ex_No_Q-(length(sg_rep_charge_ex_No_es)-1))/sg_rep_charge_ex_No_Q)*100 

sg_rep_charge_ex_Yes_es <-mods[mods$reported_charges_expunged== 'Yes', "LRIRR"]
sg_rep_charge_ex_Yes_v <-mods[mods$reported_charges_expunged== 'Yes', "RIRR_var"]
sg_rep_charge_ex_Yes_k <-length(sg_rep_charge_ex_Yes_es)
sg_rep_charge_ex_Yes_PLRIRR <-(sum(sg_rep_charge_ex_Yes_es/sg_rep_charge_ex_Yes_v))/(sum(1/sg_rep_charge_ex_Yes_v))
sg_rep_charge_ex_Yes_PLRIRR_v <- 1/sum(1/sg_rep_charge_ex_Yes_v)
sg_rep_charge_ex_Yes_Q <-sum((sg_rep_charge_ex_Yes_es-sg_rep_charge_ex_Yes_PLRIRR)^2/sg_rep_charge_ex_Yes_v)
sg_rep_charge_ex_Yes_Q_p <-pchisq(sg_rep_charge_ex_Yes_Q, df=(length(sg_rep_charge_ex_Yes_es)-1), lower.tail = FALSE)
sg_rep_charge_ex_Yes_I2 <-((sg_rep_charge_ex_Yes_Q-(length(sg_rep_charge_ex_Yes_es)-1))/sg_rep_charge_ex_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.char.ex<-lm(LRIRR ~ 1 + mods$reported_charges_expunged,weights = 1/v)

###Extract moderator Beta
mod_rep_charge_ex_LRIRR<-summary(model.lm.mod.rep.char.ex)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_charge_ex_se<-(summary(model.lm.mod.rep.char.ex)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_charge_ex_res_df<- model.lm.mod.rep.char.ex$df.residual

###Calculate Qres for the model
mod_rep_charge_ex_Y <- as.matrix(mods$LRIRR)
mod_rep_charge_ex_W <- diag(1/mods$RIRR_var)
mod_rep_charge_ex_X <- model.matrix(model.lm.mod.rep.char.ex)
mod_rep_charge_ex_P <- mod_rep_charge_ex_W - mod_rep_charge_ex_W %*% mod_rep_charge_ex_X %*% solve(
  t(mod_rep_charge_ex_X) %*% mod_rep_charge_ex_W %*% mod_rep_charge_ex_X) %*% t(
    mod_rep_charge_ex_X) %*% mod_rep_charge_ex_W 

mod_rep_charge_ex_Qres <- max(0, c(crossprod(mod_rep_charge_ex_Y, mod_rep_charge_ex_P) %*% mod_rep_charge_ex_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_charge_ex_list<-c(mod_rep_charge_ex_names,sg_rep_charge_ex_No_k, sg_rep_charge_ex_No_PLRIRR,sg_rep_charge_ex_No_PLRIRR_v, sg_rep_charge_ex_No_Q, sg_rep_charge_ex_No_Q_p, 
                          sg_rep_charge_ex_No_I2,sg_rep_charge_ex_Yes_k, sg_rep_charge_ex_Yes_PLRIRR,sg_rep_charge_ex_Yes_PLRIRR_v, sg_rep_charge_ex_Yes_Q, 
                          sg_rep_charge_ex_Yes_Q_p, sg_rep_charge_ex_Yes_I2, mod_rep_charge_ex_LRIRR, mod_rep_charge_ex_se,
                          mod_rep_charge_ex_res_df, mod_rep_charge_ex_Qres)

names(mod_rep_charge_ex_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                 "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                 "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_charge_ex_list)
########################################
### 8.4.9) Stage of disposition
########################################
###Create a new stage of disposition variable for the model
mods$stage_of_disposition_recode<- mods$stage_of_disposition
if (class(mods$stage_of_disposition_recode) == 'factor') {relevel(mods$stage_of_disposition_recode, ref = "Pre & Post")} # re-code reference level if need be

### Create an object for the row name, level_ref, and level_alternative names
mod_stag_disp_names<- c('stage_of_disposition','programmatic', 'Pre & Post', 'Post') 

###Recode Pre' as 'Pre & Post' 
mods$stage_of_disposition_recode <-recode(mods$stage_of_disposition_recode,`Pre`='Pre & Post')

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_stag_disp_pre_post_es <-mods[mods$stage_of_disposition_recode== 'Pre & Post', "LRIRR"]
sg_stag_disp_pre_post_v <-mods[mods$stage_of_disposition_recode== 'Pre & Post', "RIRR_var"]
sg_stag_disp_pre_post_k <-length(sg_stag_disp_pre_post_es)
sg_stag_disp_pre_post_PLRIRR <-(sum(sg_stag_disp_pre_post_es/sg_stag_disp_pre_post_v))/(sum(1/sg_stag_disp_pre_post_v))
sg_stag_disp_pre_post_PLRIRR_v <- 1/sum(1/sg_stag_disp_pre_post_v)
sg_stag_disp_pre_post_Q <-sum((sg_stag_disp_pre_post_es-sg_stag_disp_pre_post_PLRIRR)^2/sg_stag_disp_pre_post_v)
sg_stag_disp_pre_post_Q_p <-pchisq(sg_stag_disp_pre_post_Q, df=(length(sg_stag_disp_pre_post_es)-1), lower.tail = FALSE)
sg_stag_disp_pre_post_I2 <-((sg_stag_disp_pre_post_Q-(length(sg_stag_disp_pre_post_es)-1))/sg_stag_disp_pre_post_Q)*100

sg_stag_disp_post_only_es <-mods[mods$stage_of_disposition_recode== 'Post', "LRIRR"]
sg_stag_disp_post_only_v <-mods[mods$stage_of_disposition_recode== 'Post', "RIRR_var"]
sg_stag_disp_post_only_k <-length(sg_stag_disp_post_only_es)
sg_stag_disp_post_only_PLRIRR <-(sum(sg_stag_disp_post_only_es/sg_stag_disp_post_only_v))/(sum(1/sg_stag_disp_post_only_v))
sg_stag_disp_post_only_PLRIRR_v <- 1/sum(1/sg_stag_disp_post_only_v)
sg_stag_disp_post_only_Q <-sum((sg_stag_disp_post_only_es-sg_stag_disp_post_only_PLRIRR)^2/sg_stag_disp_post_only_v)
sg_stag_disp_post_only_Q_p <-pchisq(sg_stag_disp_post_only_Q, df=(length(sg_stag_disp_post_only_es)-1), lower.tail = FALSE)
sg_stag_disp_post_only_I2 <-((sg_stag_disp_post_only_Q-(length(sg_stag_disp_post_only_es)-1))/sg_stag_disp_post_only_Q)*100

###Recode 'Not Reported' as NA
mods$stage_of_disposition_recode <-na_if(mods$stage_of_disposition_recode, 'Not Reported')
if (class(mods$stage_of_disposition_recode) == 'factor') {droplevels(mods$stage_of_disposition_recode)}

###Insert the moderator into the model
model.lm.mod.stag.disp<-lm(LRIRR ~ 1 + mods$stage_of_disposition_recode,weights = 1/v)

###Extract moderator Beta
mod_stag_disp_LRIRR<-summary(model.lm.mod.stag.disp)$coefficients[2, 1]

###Extract the Beta's standard error
mod_stag_disp_se<-summary(model.lm.mod.stag.disp)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_stag_disp_res_df<- model.lm.mod.stag.disp$df.residual

###Calculate Qres for the model
mod_stag_disp_es <-mods[mods$stage_of_disposition_recode %in% c("Pre & Post", "Post"), "LRIRR"]
mod_stag_disp_v <-mods[mods$stage_of_disposition_recode %in% c("Pre & Post", "Post"), "RIRR_var"]

mod_stag_disp_Y <- as.matrix(mod_stag_disp_es)
mod_stag_disp_W <- diag(1/mod_stag_disp_v)
mod_stag_disp_X <- model.matrix(model.lm.mod.stag.disp)
mod_stag_disp_P <- mod_stag_disp_W - mod_stag_disp_W %*% mod_stag_disp_X %*% solve(
  t(mod_stag_disp_X) %*% mod_stag_disp_W %*% mod_stag_disp_X) %*% t(
    mod_stag_disp_X) %*% mod_stag_disp_W 

mod_stag_disp_Qres <- max(0, c(crossprod(mod_stag_disp_Y, mod_stag_disp_P) %*% mod_stag_disp_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_stag_disp_list<- c(mod_stag_disp_names,sg_stag_disp_pre_post_k,  sg_stag_disp_pre_post_PLRIRR,sg_stag_disp_pre_post_PLRIRR_v,
                       sg_stag_disp_pre_post_Q, sg_stag_disp_pre_post_Q_p, sg_stag_disp_pre_post_I2, sg_stag_disp_post_only_k, 
                       sg_stag_disp_post_only_PLRIRR,sg_stag_disp_post_only_PLRIRR_v, sg_stag_disp_post_only_Q, sg_stag_disp_post_only_Q_p,
                       sg_stag_disp_post_only_I2, mod_stag_disp_LRIRR, mod_stag_disp_se, 
                       mod_stag_disp_res_df, mod_stag_disp_Qres)

names(mod_stag_disp_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                             "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                             "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_stag_disp_list)

########################################
### 8.4.10) Single vs multiple judges presiding over court hearings
########################################
###Create a new JOs presiding variable for the model
mods$jos_presiding_recode<- mods$JOs_presiding_over_hearings_collapsed

### Create an object for the row name, level_ref, and level_alternative names
mod_jos_pres_names<- c('Judicial_officers_presiding_over_court','programmatic', 'Multiple', 'Single') 

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_jos_pres_multi_es <-mods[mods$jos_presiding_recode== 'Multiple', "LRIRR"]
sg_jos_pres_multi_v <-mods[mods$jos_presiding_recode== 'Multiple', "RIRR_var"]
sg_jos_pres_multi_k <-length(sg_jos_pres_multi_es)
sg_jos_pres_multi_PLRIRR <-(sum(sg_jos_pres_multi_es/sg_jos_pres_multi_v))/(sum(1/sg_jos_pres_multi_v))
sg_jos_pres_multi_PLRIRR_v <- 1/sum(1/sg_jos_pres_multi_v)
sg_jos_pres_multi_Q <-sum((sg_jos_pres_multi_es-sg_jos_pres_multi_PLRIRR)^2/sg_jos_pres_multi_v)
sg_jos_pres_multi_Q_p <-pchisq(sg_jos_pres_multi_Q, df=(length(sg_jos_pres_multi_es)-1), lower.tail = FALSE)
sg_jos_pres_multi_I2 <-((sg_jos_pres_multi_Q-(length(sg_jos_pres_multi_es)-1))/sg_jos_pres_multi_Q)*100

sg_jos_pres_single_es <-mods[mods$jos_presiding_recode== 'Single judge', "LRIRR"]
sg_jos_pres_single_v <-mods[mods$jos_presiding_recode== 'Single judge', "RIRR_var"]
sg_jos_pres_single_k <-length(sg_jos_pres_single_es)
sg_jos_pres_single_PLRIRR <-(sum(sg_jos_pres_single_es/sg_jos_pres_single_v))/(sum(1/sg_jos_pres_single_v))
sg_jos_pres_single_PLRIRR_v <- 1/sum(1/sg_jos_pres_single_v)
sg_jos_pres_single_Q <-sum((sg_jos_pres_single_es-sg_jos_pres_single_PLRIRR)^2/sg_jos_pres_single_v)
sg_jos_pres_single_Q_p <-pchisq(sg_jos_pres_single_Q, df=(length(sg_jos_pres_single_es)-1), lower.tail = FALSE)
sg_jos_pres_single_I2 <-((sg_jos_pres_single_Q-(length(sg_jos_pres_single_es)-1))/sg_jos_pres_single_Q)*100

###Recode 'Not Reported' as NA
mods$jos_presiding_recode <-na_if(mods$jos_presiding_recode, 'Not Reported')
# mods$jos_presiding_recode <-droplevels(mods$jos_presiding_recode)

###Insert the moderator into the model
model.lm.mod.jos.pres<-lm(LRIRR ~ 1 + mods$jos_presiding_recode,weights = 1/v)

###Extract moderator Beta
mod_jos_pres_LRIRR<-summary(model.lm.mod.jos.pres)$coefficients[2, 1]

###Extract the Beta's standard error
mod_jos_pres_se<-summary(model.lm.mod.jos.pres)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_jos_pres_res_df<- model.lm.mod.jos.pres$df.residual

###Calculate Qres for the model
mod_jos_pres_es <-mods[mods$jos_presiding_recode %in% c("Multiple", "Single judge"), "LRIRR"]
mod_jos_pres_v <-mods[mods$jos_presiding_recode %in% c("Multiple", "Single judge"), "RIRR_var"]

mod_jos_pres_Y <- as.matrix(mod_jos_pres_es)
mod_jos_pres_W <- diag(1/mod_jos_pres_v)
mod_jos_pres_X <- model.matrix(model.lm.mod.jos.pres)
mod_jos_pres_P <- mod_jos_pres_W - mod_jos_pres_W %*% mod_jos_pres_X %*% solve(
  t(mod_jos_pres_X) %*% mod_jos_pres_W %*% mod_jos_pres_X) %*% t(
    mod_jos_pres_X) %*% mod_jos_pres_W 

mod_jos_pres_Qres <- max(0, c(crossprod(mod_jos_pres_Y, mod_jos_pres_P) %*% mod_jos_pres_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_jos_pres_list<- c(mod_jos_pres_names,sg_jos_pres_multi_k,  sg_jos_pres_multi_PLRIRR,sg_jos_pres_multi_PLRIRR_v, 
                      sg_jos_pres_multi_Q, sg_jos_pres_multi_Q_p, sg_jos_pres_multi_I2, sg_jos_pres_single_k, 
                      sg_jos_pres_single_PLRIRR,sg_jos_pres_single_PLRIRR_v, sg_jos_pres_single_Q, sg_jos_pres_single_Q_p,
                      sg_jos_pres_single_I2, mod_jos_pres_LRIRR, mod_jos_pres_se, 
                      mod_jos_pres_res_df, mod_jos_pres_Qres)

names(mod_jos_pres_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                            "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                            "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_jos_pres_list)
########################################
### 8.4.11) Reported that a single judge presided over court hearings 
########################################
### Create an object for the row name, level_ref, and level_alternative names
reported_single_judge_names<- c('reported_single_judge_presiding', 'programmatic', 'No', 'Yes') 

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_rep_single_no_es <-mods[mods$reported_single_judge == 'No', "LRIRR"]
sg_rep_single_no_v <-mods[mods$reported_single_judge == 'No', "RIRR_var"]
sg_rep_single_no_k <-length(sg_rep_single_no_es)
sg_rep_single_no_PLRIRR <-(sum(sg_rep_single_no_es/sg_rep_single_no_v))/(sum(1/sg_rep_single_no_v))
sg_rep_single_no_PLRIRR_v <- 1/sum(1/sg_rep_single_no_v)
sg_rep_single_no_Q <-sum((sg_rep_single_no_es-sg_rep_single_no_PLRIRR)^2/sg_rep_single_no_v)
sg_rep_single_no_Q_p <-pchisq(sg_rep_single_no_Q, df=(length(sg_rep_single_no_es)-1), lower.tail = FALSE)
sg_rep_single_no_I2 <-((sg_rep_single_no_Q-(length(sg_rep_single_no_es)-1))/sg_rep_single_no_Q)*100  

sg_rep_single_yes_es <-mods[mods$reported_single_judge == 'Yes', "LRIRR"]
sg_rep_single_yes_v <-mods[mods$reported_single_judge == 'Yes', "RIRR_var"]
sg_rep_single_yes_k <-length(sg_rep_single_yes_es)
sg_rep_single_yes_PLRIRR <-sum(sg_rep_single_yes_es/sg_rep_single_yes_v)/(sum(1/sg_rep_single_yes_v))
sg_rep_single_yes_PLRIRR_v <- 1/sum(1/sg_rep_single_yes_v)
sg_rep_single_yes_Q <-sum((sg_rep_single_yes_es-sg_rep_single_yes_PLRIRR)^2/sg_rep_single_yes_v)
sg_rep_single_yes_Q_p <-pchisq(sg_rep_single_yes_Q, df=(length(sg_rep_single_yes_es)-1), lower.tail = FALSE)
sg_rep_single_yes_I2 <-((sg_rep_single_yes_Q-(length(sg_rep_single_yes_es)-1))/sg_rep_single_yes_Q)*100  

###Insert the moderator into the model
model.lm.mod.rep.single<-lm(LRIRR ~ 1 + mods$reported_single_judge,weights = 1/v)
summary(model.lm.mod.rep.single)

###Extract moderator Beta
mod_rep_single_B_PLRIRR<-summary(model.lm.mod.rep.single)$coefficients[2, 1]
mod_rep_single_B_PLRIRR

###Extract the Beta's standard error
mod_rep_single_B_PLRIRR_se<-(summary(model.lm.mod.rep.single)$coefficients[2, 2])
mod_rep_single_B_PLRIRR_se

###Create an object for the residual degrees of freedom. 
mod_rep_single_res_df<- model.lm.mod.rep.single$df.residual
mod_rep_single_res_df

###Calculate Qres for the model
###The following code is derived from Viectbauer's (2010) 'rma.uni' source code and calculates Qres for the model
mod_rep_single_Y <- as.matrix(mods$LRIRR)
mod_rep_single_W <- diag(1/v)
mod_rep_single_X <- model.matrix(model.lm.mod.rep.single)
mod_rep_single_P <- mod_rep_single_W - mod_rep_single_W %*% mod_rep_single_X %*% solve(
  t(mod_rep_single_X) %*% mod_rep_single_W %*% mod_rep_single_X) %*% t(
    mod_rep_single_X) %*% mod_rep_single_W 

mod_rep_single_Qres <- max(0, c(crossprod(mod_rep_single_Y, mod_rep_single_P) %*% mod_rep_single_Y))

###Bind the list of values to the dataframe 'moderator_output'
reported_single_judge<- c(reported_single_judge_names,sg_rep_single_no_k, sg_rep_single_no_PLRIRR,sg_rep_single_no_PLRIRR_v,
                          sg_rep_single_no_Q, sg_rep_single_no_Q_p, sg_rep_single_no_I2, sg_rep_single_yes_k,
                          sg_rep_single_yes_PLRIRR,sg_rep_single_yes_PLRIRR_v, sg_rep_single_yes_Q, sg_rep_single_yes_Q_p, 
                          sg_rep_single_yes_I2, mod_rep_single_B_PLRIRR, mod_rep_single_B_PLRIRR_se, 
                          mod_rep_single_res_df, mod_rep_single_Qres)

names(reported_single_judge)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                "mod_Qres")
moderator_output<- bind_rows(moderator_output, reported_single_judge)
########################################
### 8.4.12) Reported team training 
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_team_train_names<- c('reported_team_training','programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_team_train_No_es <-mods[mods$report_team_training== 'No', "LRIRR"]
sg_rep_team_train_No_v <-mods[mods$report_team_training== 'No', "RIRR_var"]
sg_rep_team_train_No_k <-length(sg_rep_team_train_No_es)
sg_rep_team_train_No_PLRIRR <-(sum(sg_rep_team_train_No_es/sg_rep_team_train_No_v))/(sum(1/sg_rep_team_train_No_v))
sg_rep_team_train_No_PLRIRR_v <- 1/sum(1/sg_rep_team_train_No_v)
sg_rep_team_train_No_Q <-sum((sg_rep_team_train_No_es-sg_rep_team_train_No_PLRIRR)^2/sg_rep_team_train_No_v)
sg_rep_team_train_No_Q_p <-pchisq(sg_rep_team_train_No_Q, df=(length(sg_rep_team_train_No_es)-1), lower.tail = FALSE)
sg_rep_team_train_No_I2 <-((sg_rep_team_train_No_Q-(length(sg_rep_team_train_No_es)-1))/sg_rep_team_train_No_Q)*100 

sg_rep_team_train_Yes_es <-mods[mods$report_team_training== 'Yes', "LRIRR"]
sg_rep_team_train_Yes_v <-mods[mods$report_team_training== 'Yes', "RIRR_var"]
sg_rep_team_train_Yes_k <-length(sg_rep_team_train_Yes_es)
sg_rep_team_train_Yes_PLRIRR <-(sum(sg_rep_team_train_Yes_es/sg_rep_team_train_Yes_v))/(sum(1/sg_rep_team_train_Yes_v))
sg_rep_team_train_Yes_PLRIRR_v <- 1/sum(1/sg_rep_team_train_Yes_v)
sg_rep_team_train_Yes_Q <-sum((sg_rep_team_train_Yes_es-sg_rep_team_train_Yes_PLRIRR)^2/sg_rep_team_train_Yes_v)
sg_rep_team_train_Yes_Q_p <-pchisq(sg_rep_team_train_Yes_Q, df=(length(sg_rep_team_train_Yes_es)-1), lower.tail = FALSE)
sg_rep_team_train_Yes_I2 <-((sg_rep_team_train_Yes_Q-(length(sg_rep_team_train_Yes_es)-1))/sg_rep_team_train_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.team.train<-lm(LRIRR ~ 1 + mods$report_team_training,weights = 1/v)

###Extract moderator Beta
mod_rep_team_train_LRIRR<-summary(model.lm.mod.rep.team.train)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_team_train_se<-(summary(model.lm.mod.rep.team.train)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_team_train_res_df<- model.lm.mod.rep.team.train$df.residual

###Calculate Qres for the model
mod_rep_team_train_Y <- as.matrix(mods$LRIRR)
mod_rep_team_train_W <- diag(1/mods$RIRR_var)
mod_rep_team_train_X <- model.matrix(model.lm.mod.rep.team.train)
mod_rep_team_train_P <- mod_rep_team_train_W - mod_rep_team_train_W %*% mod_rep_team_train_X %*% solve(
  t(mod_rep_team_train_X) %*% mod_rep_team_train_W %*% mod_rep_team_train_X) %*% t(
    mod_rep_team_train_X) %*% mod_rep_team_train_W 

mod_rep_team_train_Qres <- max(0, c(crossprod(mod_rep_team_train_Y, mod_rep_team_train_P) %*% mod_rep_team_train_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_team_train_list<-c(mod_rep_team_train_names,sg_rep_team_train_No_k, sg_rep_team_train_No_PLRIRR,sg_rep_team_train_No_PLRIRR_v, sg_rep_team_train_No_Q, sg_rep_team_train_No_Q_p, 
                           sg_rep_team_train_No_I2,sg_rep_team_train_Yes_k, sg_rep_team_train_Yes_PLRIRR,sg_rep_team_train_Yes_PLRIRR_v, sg_rep_team_train_Yes_Q, 
                           sg_rep_team_train_Yes_Q_p, sg_rep_team_train_Yes_I2, mod_rep_team_train_LRIRR, mod_rep_team_train_se,
                           mod_rep_team_train_res_df, mod_rep_team_train_Qres)

names(mod_rep_team_train_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                  "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                  "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_team_train_list)
########################################
### 8.4.13) Reported training judge
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_jo_train_names<- c('reported_JO_training','programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_jo_train_No_es <-mods[mods$report_JO.training== 'No', "LRIRR"]
sg_rep_jo_train_No_v <-mods[mods$report_JO.training== 'No', "RIRR_var"]
sg_rep_jo_train_No_k <-length(sg_rep_jo_train_No_es)
sg_rep_jo_train_No_PLRIRR <-(sum(sg_rep_jo_train_No_es/sg_rep_jo_train_No_v))/(sum(1/sg_rep_jo_train_No_v))
sg_rep_jo_train_No_PLRIRR_v <- 1/sum(1/sg_rep_jo_train_No_v)
sg_rep_jo_train_No_Q <-sum((sg_rep_jo_train_No_es-sg_rep_jo_train_No_PLRIRR)^2/sg_rep_jo_train_No_v)
sg_rep_jo_train_No_Q_p <-pchisq(sg_rep_jo_train_No_Q, df=(length(sg_rep_jo_train_No_es)-1), lower.tail = FALSE)
sg_rep_jo_train_No_I2 <-((sg_rep_jo_train_No_Q-(length(sg_rep_jo_train_No_es)-1))/sg_rep_jo_train_No_Q)*100 

sg_rep_jo_train_Yes_es <-mods[mods$report_JO.training== 'Yes', "LRIRR"]
sg_rep_jo_train_Yes_v <-mods[mods$report_JO.training== 'Yes', "RIRR_var"]
sg_rep_jo_train_Yes_k <-length(sg_rep_jo_train_Yes_es)
sg_rep_jo_train_Yes_PLRIRR <-(sum(sg_rep_jo_train_Yes_es/sg_rep_jo_train_Yes_v))/(sum(1/sg_rep_jo_train_Yes_v))
sg_rep_jo_train_Yes_PLRIRR_v <- 1/sum(1/sg_rep_jo_train_Yes_v)
sg_rep_jo_train_Yes_Q <-sum((sg_rep_jo_train_Yes_es-sg_rep_jo_train_Yes_PLRIRR)^2/sg_rep_jo_train_Yes_v)
sg_rep_jo_train_Yes_Q_p <-pchisq(sg_rep_jo_train_Yes_Q, df=(length(sg_rep_jo_train_Yes_es)-1), lower.tail = FALSE)
sg_rep_jo_train_Yes_I2 <-((sg_rep_jo_train_Yes_Q-(length(sg_rep_jo_train_Yes_es)-1))/sg_rep_jo_train_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.jo.train<-lm(LRIRR ~ 1 + mods$report_JO.training,weights = 1/v)

###Extract moderator Beta
mod_rep_jo_train_LRIRR<-summary(model.lm.mod.rep.jo.train)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_jo_train_se<-(summary(model.lm.mod.rep.jo.train)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_jo_train_res_df<- model.lm.mod.rep.jo.train$df.residual

###Calculate Qres for the model
mod_rep_jo_train_Y <- as.matrix(mods$LRIRR)
mod_rep_jo_train_W <- diag(1/mods$RIRR_var)
mod_rep_jo_train_X <- model.matrix(model.lm.mod.rep.jo.train)
mod_rep_jo_train_P <- mod_rep_jo_train_W - mod_rep_jo_train_W %*% mod_rep_jo_train_X %*% solve(
  t(mod_rep_jo_train_X) %*% mod_rep_jo_train_W %*% mod_rep_jo_train_X) %*% t(
    mod_rep_jo_train_X) %*% mod_rep_jo_train_W 

mod_rep_jo_train_Qres <- max(0, c(crossprod(mod_rep_jo_train_Y, mod_rep_jo_train_P) %*% mod_rep_jo_train_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_jo_train_list<-c(mod_rep_jo_train_names,sg_rep_jo_train_No_k, sg_rep_jo_train_No_PLRIRR,sg_rep_jo_train_No_PLRIRR_v, sg_rep_jo_train_No_Q, sg_rep_jo_train_No_Q_p, 
                         sg_rep_jo_train_No_I2,sg_rep_jo_train_Yes_k, sg_rep_jo_train_Yes_PLRIRR,sg_rep_jo_train_Yes_PLRIRR_v, sg_rep_jo_train_Yes_Q, 
                         sg_rep_jo_train_Yes_Q_p, sg_rep_jo_train_Yes_I2, mod_rep_jo_train_LRIRR, mod_rep_jo_train_se,
                         mod_rep_jo_train_res_df, mod_rep_jo_train_Qres)

names(mod_rep_jo_train_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_jo_train_list)
########################################
### 8.4.14) Reported training judge or team
########################################
# Note: this variable separates the subgroups in the same way as '8.4.12) Reported training team' and so will be excluded from the mean moderator analysis.
all.equal(table(mods$report_judge_or_team_training, mods$studyid), table(mods$report_team_training, mods$studyid))
########################################
### 8.4.15) Average minutes per court review hearing reported for court
########################################
avg_min_names<- c('average_review_hearing(minutes)', 'programmatic', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                  'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                  'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.avg.min<-lm(LRIRR ~ 1 + mods$Average._minutes_review_hearing, weights = 1/v)

###Extract moderator Beta
mod_avg_min_B_LRIRR<-summary(model.lm.mod.avg.min)$coefficients[2, 1]

###Extract the Beta's standard error
mod_avg_min_B_se<-summary(model.lm.mod.avg.min)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_avg_min_res_df<- model.lm.mod.avg.min$df.residual

###Calculate Qres for the model
mod_mod_avg_es <-mods[mods$Average._minutes_review_hearing %in% na.omit(mods$Average._minutes_review_hearing), "LRIRR"]
mod_mod_avg_v <-mods[mods$Average._minutes_review_hearing %in% na.omit(mods$Average._minutes_review_hearing), "RIRR_var"]

mod_avg_min_Y <- as.matrix(mod_mod_avg_es)
mod_avg_min_W <- diag(1/mod_mod_avg_v)
mod_avg_min_X <- model.matrix(model.lm.mod.avg.min)
mod_avg_min_P <- mod_avg_min_W - mod_avg_min_W %*% mod_avg_min_X %*% solve(
  t(mod_avg_min_X) %*% mod_avg_min_W %*% mod_avg_min_X) %*% t(
    mod_avg_min_X) %*% mod_avg_min_W 

mod_avg_min_Qres <- max(0, c(crossprod(mod_avg_min_Y, mod_avg_min_P) %*% mod_avg_min_Y))

###Bind the list of values to the dataframe 'moderator_output'
avg_min_name_list<- c(avg_min_names, mod_avg_min_B_LRIRR, mod_avg_min_B_se, 
                      mod_avg_min_res_df, mod_avg_min_Qres)

names(avg_min_name_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                            "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                            "mod_Qres")
moderator_output<- bind_rows(moderator_output, avg_min_name_list)
########################################
### 8.4.16) Reported frequent team meetings
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_team_meet_names<- c('reported_frequent_team_meetings','programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_team_meet_No_es <-mods[mods$report_frequent_team_meetings== 'No', "LRIRR"]
sg_rep_team_meet_No_v <-mods[mods$report_frequent_team_meetings== 'No', "RIRR_var"]
sg_rep_team_meet_No_k <-length(sg_rep_team_meet_No_es)
sg_rep_team_meet_No_PLRIRR <-(sum(sg_rep_team_meet_No_es/sg_rep_team_meet_No_v))/(sum(1/sg_rep_team_meet_No_v))
sg_rep_team_meet_No_PLRIRR_v <- 1/sum(1/sg_rep_team_meet_No_v)
sg_rep_team_meet_No_Q <-sum((sg_rep_team_meet_No_es-sg_rep_team_meet_No_PLRIRR)^2/sg_rep_team_meet_No_v)
sg_rep_team_meet_No_Q_p <-pchisq(sg_rep_team_meet_No_Q, df=(length(sg_rep_team_meet_No_es)-1), lower.tail = FALSE)
sg_rep_team_meet_No_I2 <-((sg_rep_team_meet_No_Q-(length(sg_rep_team_meet_No_es)-1))/sg_rep_team_meet_No_Q)*100 

sg_rep_team_meet_Yes_es <-mods[mods$report_frequent_team_meetings== 'Yes', "LRIRR"]
sg_rep_team_meet_Yes_v <-mods[mods$report_frequent_team_meetings== 'Yes', "RIRR_var"]
sg_rep_team_meet_Yes_k <-length(sg_rep_team_meet_Yes_es)
sg_rep_team_meet_Yes_PLRIRR <-(sum(sg_rep_team_meet_Yes_es/sg_rep_team_meet_Yes_v))/(sum(1/sg_rep_team_meet_Yes_v))
sg_rep_team_meet_Yes_PLRIRR_v <- 1/sum(1/sg_rep_team_meet_Yes_v)
sg_rep_team_meet_Yes_Q <-sum((sg_rep_team_meet_Yes_es-sg_rep_team_meet_Yes_PLRIRR)^2/sg_rep_team_meet_Yes_v)
sg_rep_team_meet_Yes_Q_p <-pchisq(sg_rep_team_meet_Yes_Q, df=(length(sg_rep_team_meet_Yes_es)-1), lower.tail = FALSE)
sg_rep_team_meet_Yes_I2 <-((sg_rep_team_meet_Yes_Q-(length(sg_rep_team_meet_Yes_es)-1))/sg_rep_team_meet_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.team.meet<-lm(LRIRR ~ 1 + mods$report_frequent_team_meetings,weights = 1/v)

###Extract moderator Beta
mod_rep_team_meet_LRIRR<-summary(model.lm.mod.rep.team.meet)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_team_meet_se<-(summary(model.lm.mod.rep.team.meet)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_team_meet_res_df<- model.lm.mod.rep.team.meet$df.residual

###Calculate Qres for the model
mod_rep_team_meet_Y <- as.matrix(mods$LRIRR)
mod_rep_team_meet_W <- diag(1/mods$RIRR_var)
mod_rep_team_meet_X <- model.matrix(model.lm.mod.rep.team.meet)
mod_rep_team_meet_P <- mod_rep_team_meet_W - mod_rep_team_meet_W %*% mod_rep_team_meet_X %*% solve(
  t(mod_rep_team_meet_X) %*% mod_rep_team_meet_W %*% mod_rep_team_meet_X) %*% t(
    mod_rep_team_meet_X) %*% mod_rep_team_meet_W 

mod_rep_team_meet_Qres <- max(0, c(crossprod(mod_rep_team_meet_Y, mod_rep_team_meet_P) %*% mod_rep_team_meet_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_team_meet_list<-c(mod_rep_team_meet_names,sg_rep_team_meet_No_k, sg_rep_team_meet_No_PLRIRR,sg_rep_team_meet_No_PLRIRR_v, sg_rep_team_meet_No_Q, sg_rep_team_meet_No_Q_p, 
                          sg_rep_team_meet_No_I2,sg_rep_team_meet_Yes_k, sg_rep_team_meet_Yes_PLRIRR,sg_rep_team_meet_Yes_PLRIRR_v, sg_rep_team_meet_Yes_Q, 
                          sg_rep_team_meet_Yes_Q_p, sg_rep_team_meet_Yes_I2, mod_rep_team_meet_LRIRR, mod_rep_team_meet_se,
                          mod_rep_team_meet_res_df, mod_rep_team_meet_Qres)

names(mod_rep_team_meet_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                 "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                 "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_team_meet_list)
########################################
### 8.4.17) Reported using sanctions
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_sanctions_names<- c('reported_sanctions','programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_sanctions_No_es <-mods[mods$reported_sanctions== 'No', "LRIRR"]
sg_rep_sanctions_No_v <-mods[mods$reported_sanctions== 'No', "RIRR_var"]
sg_rep_sanctions_No_k <-length(sg_rep_sanctions_No_es)
sg_rep_sanctions_No_PLRIRR <-(sum(sg_rep_sanctions_No_es/sg_rep_sanctions_No_v))/(sum(1/sg_rep_sanctions_No_v))
sg_rep_sanctions_No_PLRIRR_v <- 1/sum(1/sg_rep_sanctions_No_v)
sg_rep_sanctions_No_Q <-sum((sg_rep_sanctions_No_es-sg_rep_sanctions_No_PLRIRR)^2/sg_rep_sanctions_No_v)
sg_rep_sanctions_No_Q_p <-pchisq(sg_rep_sanctions_No_Q, df=(length(sg_rep_sanctions_No_es)-1), lower.tail = FALSE)
sg_rep_sanctions_No_I2 <-((sg_rep_sanctions_No_Q-(length(sg_rep_sanctions_No_es)-1))/sg_rep_sanctions_No_Q)*100 

sg_rep_sanctions_Yes_es <-mods[mods$reported_sanctions== 'Yes', "LRIRR"]
sg_rep_sanctions_Yes_v <-mods[mods$reported_sanctions== 'Yes', "RIRR_var"]
sg_rep_sanctions_Yes_k <-length(sg_rep_sanctions_Yes_es)
sg_rep_sanctions_Yes_PLRIRR <-(sum(sg_rep_sanctions_Yes_es/sg_rep_sanctions_Yes_v))/(sum(1/sg_rep_sanctions_Yes_v))
sg_rep_sanctions_Yes_PLRIRR_v <- 1/sum(1/sg_rep_sanctions_Yes_v)
sg_rep_sanctions_Yes_Q <-sum((sg_rep_sanctions_Yes_es-sg_rep_sanctions_Yes_PLRIRR)^2/sg_rep_sanctions_Yes_v)
sg_rep_sanctions_Yes_Q_p <-pchisq(sg_rep_sanctions_Yes_Q, df=(length(sg_rep_sanctions_Yes_es)-1), lower.tail = FALSE)
sg_rep_sanctions_Yes_I2 <-((sg_rep_sanctions_Yes_Q-(length(sg_rep_sanctions_Yes_es)-1))/sg_rep_sanctions_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.sanctions<-lm(LRIRR ~ 1 + mods$reported_sanctions,weights = 1/v)

###Extract moderator Beta
mod_rep_sanctions_LRIRR<-summary(model.lm.mod.rep.sanctions)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_sanctions_se<-(summary(model.lm.mod.rep.sanctions)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_sanctions_res_df<- model.lm.mod.rep.sanctions$df.residual

###Calculate Qres for the model
mod_rep_sanctions_Y <- as.matrix(mods$LRIRR)
mod_rep_sanctions_W <- diag(1/mods$RIRR_var)
mod_rep_sanctions_X <- model.matrix(model.lm.mod.rep.sanctions)
mod_rep_sanctions_P <- mod_rep_sanctions_W - mod_rep_sanctions_W %*% mod_rep_sanctions_X %*% solve(
  t(mod_rep_sanctions_X) %*% mod_rep_sanctions_W %*% mod_rep_sanctions_X) %*% t(
    mod_rep_sanctions_X) %*% mod_rep_sanctions_W 

mod_rep_sanctions_Qres <- max(0, c(crossprod(mod_rep_sanctions_Y, mod_rep_sanctions_P) %*% mod_rep_sanctions_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_sanctions_list<-c(mod_rep_sanctions_names,sg_rep_sanctions_No_k, sg_rep_sanctions_No_PLRIRR,sg_rep_sanctions_No_PLRIRR_v, sg_rep_sanctions_No_Q, sg_rep_sanctions_No_Q_p, 
                          sg_rep_sanctions_No_I2,sg_rep_sanctions_Yes_k, sg_rep_sanctions_Yes_PLRIRR,sg_rep_sanctions_Yes_PLRIRR_v, sg_rep_sanctions_Yes_Q, 
                          sg_rep_sanctions_Yes_Q_p, sg_rep_sanctions_Yes_I2, mod_rep_sanctions_LRIRR, mod_rep_sanctions_se,
                          mod_rep_sanctions_res_df, mod_rep_sanctions_Qres)

names(mod_rep_sanctions_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                 "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                 "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_sanctions_list)
########################################
### 8.4.18) Reported using graduated sanctions
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_grad_sanctions_names<- c('reported_graduated_sanctions','programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_grad_sanctions_No_es <-mods[mods$reported_using_graduated_sanctions== 'No', "LRIRR"]
sg_rep_grad_sanctions_No_v <-mods[mods$reported_using_graduated_sanctions== 'No', "RIRR_var"]
sg_rep_grad_sanctions_No_k <-length(sg_rep_grad_sanctions_No_es)
sg_rep_grad_sanctions_No_PLRIRR <-(sum(sg_rep_grad_sanctions_No_es/sg_rep_grad_sanctions_No_v))/(sum(1/sg_rep_grad_sanctions_No_v))
sg_rep_grad_sanctions_No_PLRIRR_v <- 1/sum(1/sg_rep_grad_sanctions_No_v)
sg_rep_grad_sanctions_No_Q <-sum((sg_rep_grad_sanctions_No_es-sg_rep_grad_sanctions_No_PLRIRR)^2/sg_rep_grad_sanctions_No_v)
sg_rep_grad_sanctions_No_Q_p <-pchisq(sg_rep_grad_sanctions_No_Q, df=(length(sg_rep_grad_sanctions_No_es)-1), lower.tail = FALSE)
sg_rep_grad_sanctions_No_I2 <-((sg_rep_grad_sanctions_No_Q-(length(sg_rep_grad_sanctions_No_es)-1))/sg_rep_grad_sanctions_No_Q)*100 

sg_rep_grad_sanctions_Yes_es <-mods[mods$reported_using_graduated_sanctions== 'Yes', "LRIRR"]
sg_rep_grad_sanctions_Yes_v <-mods[mods$reported_using_graduated_sanctions== 'Yes', "RIRR_var"]
sg_rep_grad_sanctions_Yes_k <-length(sg_rep_grad_sanctions_Yes_es)
sg_rep_grad_sanctions_Yes_PLRIRR <-(sum(sg_rep_grad_sanctions_Yes_es/sg_rep_grad_sanctions_Yes_v))/(sum(1/sg_rep_grad_sanctions_Yes_v))
sg_rep_grad_sanctions_Yes_PLRIRR_v <- 1/sum(1/sg_rep_grad_sanctions_Yes_v)
sg_rep_grad_sanctions_Yes_Q <-sum((sg_rep_grad_sanctions_Yes_es-sg_rep_grad_sanctions_Yes_PLRIRR)^2/sg_rep_grad_sanctions_Yes_v)
sg_rep_grad_sanctions_Yes_Q_p <-pchisq(sg_rep_grad_sanctions_Yes_Q, df=(length(sg_rep_grad_sanctions_Yes_es)-1), lower.tail = FALSE)
sg_rep_grad_sanctions_Yes_I2 <-((sg_rep_grad_sanctions_Yes_Q-(length(sg_rep_grad_sanctions_Yes_es)-1))/sg_rep_grad_sanctions_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.grad.sanctions<-lm(LRIRR ~ 1 + mods$reported_using_graduated_sanctions,weights = 1/v)

###Extract moderator Beta
mod_rep_grad_sanctions_LRIRR<-summary(model.lm.mod.rep.grad.sanctions)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_grad_sanctions_se<-(summary(model.lm.mod.rep.grad.sanctions)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_grad_sanctions_res_df<- model.lm.mod.rep.grad.sanctions$df.residual

###Calculate Qres for the model
mod_rep_grad_sanctions_Y <- as.matrix(mods$LRIRR)
mod_rep_grad_sanctions_W <- diag(1/mods$RIRR_var)
mod_rep_grad_sanctions_X <- model.matrix(model.lm.mod.rep.grad.sanctions)
mod_rep_grad_sanctions_P <- mod_rep_grad_sanctions_W - mod_rep_grad_sanctions_W %*% mod_rep_grad_sanctions_X %*% solve(
  t(mod_rep_grad_sanctions_X) %*% mod_rep_grad_sanctions_W %*% mod_rep_grad_sanctions_X) %*% t(
    mod_rep_grad_sanctions_X) %*% mod_rep_grad_sanctions_W 

mod_rep_grad_sanctions_Qres <- max(0, c(crossprod(mod_rep_grad_sanctions_Y, mod_rep_grad_sanctions_P) %*% mod_rep_grad_sanctions_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_grad_sanctions_list<-c(mod_rep_grad_sanctions_names,sg_rep_grad_sanctions_No_k, sg_rep_grad_sanctions_No_PLRIRR,sg_rep_grad_sanctions_No_PLRIRR_v, sg_rep_grad_sanctions_No_Q, sg_rep_grad_sanctions_No_Q_p, 
                               sg_rep_grad_sanctions_No_I2,sg_rep_grad_sanctions_Yes_k, sg_rep_grad_sanctions_Yes_PLRIRR,sg_rep_grad_sanctions_Yes_PLRIRR_v, sg_rep_grad_sanctions_Yes_Q, 
                               sg_rep_grad_sanctions_Yes_Q_p, sg_rep_grad_sanctions_Yes_I2, mod_rep_grad_sanctions_LRIRR, mod_rep_grad_sanctions_se,
                               mod_rep_grad_sanctions_res_df, mod_rep_grad_sanctions_Qres)

names(mod_rep_grad_sanctions_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                      "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                      "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_grad_sanctions_list)
########################################
### 8.4.19.1) Reported sanctions or graduated sanctions - Graduated
########################################
###Re-code the setting variable for the model if need be
if (class(mods$reported_sanctions_or_graduated_sanctions) == 'factor') {relevel(mods$reported_sanctions_or_graduated_sanctions, ref = "No")} 

### Create an object for the row name, level_ref, and level_alternative names
mod_rep_grad_or_sanctions_grad_names<- c('reported_sanctions_or_graduated_sanctions','programmatic', 'No', 'graduated')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_grad_or_sanctions_No_es <-mods[mods$reported_sanctions_or_graduated_sanctions== 'No', "LRIRR"]
sg_rep_grad_or_sanctions_No_v <-mods[mods$reported_sanctions_or_graduated_sanctions== 'No', "RIRR_var"]
sg_rep_grad_or_sanctions_No_k <-length(sg_rep_grad_or_sanctions_No_es)
sg_rep_grad_or_sanctions_No_PLRIRR <-(sum(sg_rep_grad_or_sanctions_No_es/sg_rep_grad_or_sanctions_No_v))/(sum(1/sg_rep_grad_or_sanctions_No_v))
sg_rep_grad_or_sanctions_No_PLRIRR_v <- 1/sum(1/sg_rep_grad_or_sanctions_No_v)
sg_rep_grad_or_sanctions_No_Q <-sum((sg_rep_grad_or_sanctions_No_es-sg_rep_grad_or_sanctions_No_PLRIRR)^2/sg_rep_grad_or_sanctions_No_v)
sg_rep_grad_or_sanctions_No_Q_p <-pchisq(sg_rep_grad_or_sanctions_No_Q, df=(length(sg_rep_grad_or_sanctions_No_es)-1), lower.tail = FALSE)
sg_rep_grad_or_sanctions_No_I2 <-((sg_rep_grad_or_sanctions_No_Q-(length(sg_rep_grad_or_sanctions_No_es)-1))/sg_rep_grad_or_sanctions_No_Q)*100 

sg_rep_grad_or_sanctions_grad_es <-mods[mods$reported_sanctions_or_graduated_sanctions== 'Graduated', "LRIRR"]
sg_rep_grad_or_sanctions_grad_v <-mods[mods$reported_sanctions_or_graduated_sanctions== 'Graduated', "RIRR_var"]
sg_rep_grad_or_sanctions_grad_k <-length(sg_rep_grad_or_sanctions_grad_es)
sg_rep_grad_or_sanctions_grad_PLRIRR <-(sum(sg_rep_grad_or_sanctions_grad_es/sg_rep_grad_or_sanctions_grad_v))/(sum(1/sg_rep_grad_or_sanctions_grad_v))
sg_rep_grad_or_sanctions_grad_PLRIRR_v <- 1/sum(1/sg_rep_grad_or_sanctions_grad_v)
sg_rep_grad_or_sanctions_grad_Q <-sum((sg_rep_grad_or_sanctions_grad_es-sg_rep_grad_or_sanctions_grad_PLRIRR)^2/sg_rep_grad_or_sanctions_grad_v)
sg_rep_grad_or_sanctions_grad_Q_p <-pchisq(sg_rep_grad_or_sanctions_grad_Q, df=(length(sg_rep_grad_or_sanctions_grad_es)-1), lower.tail = FALSE)
sg_rep_grad_or_sanctions_grad_I2 <-((sg_rep_grad_or_sanctions_grad_Q-(length(sg_rep_grad_or_sanctions_grad_es)-1))/sg_rep_grad_or_sanctions_grad_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.grad.or.sanct<-lm(LRIRR ~ 1 + mods$reported_sanctions_or_graduated_sanctions,weights = 1/v)

###Extract moderator Beta
mod_rep_grad_or_sanctions_LRIRR<-summary(model.lm.mod.rep.grad.or.sanct)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_grad_or_sanctions_se<-(summary(model.lm.mod.rep.grad.or.sanct)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_grad_or_sanctions_res_df<- model.lm.mod.rep.grad.or.sanct$df.residual

###Calculate Qres for the model
mod_rep_grad_or_sanctions_Y <- as.matrix(mods$LRIRR)
mod_rep_grad_or_sanctions_W <- diag(1/mods$RIRR_var)
mod_rep_grad_or_sanctions_X <- model.matrix(model.lm.mod.rep.grad.or.sanct)
mod_rep_grad_or_sanctions_P <- mod_rep_grad_or_sanctions_W - mod_rep_grad_or_sanctions_W %*% mod_rep_grad_or_sanctions_X %*% solve(
  t(mod_rep_grad_or_sanctions_X) %*% mod_rep_grad_or_sanctions_W %*% mod_rep_grad_or_sanctions_X) %*% t(
    mod_rep_grad_or_sanctions_X) %*% mod_rep_grad_or_sanctions_W 

mod_rep_grad_or_sanctions_Qres <- max(0, c(crossprod(mod_rep_grad_or_sanctions_Y, mod_rep_grad_or_sanctions_P) %*% mod_rep_grad_or_sanctions_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_grad_or_sanctions_list<-c(mod_rep_grad_or_sanctions_grad_names,sg_rep_grad_or_sanctions_No_k, sg_rep_grad_or_sanctions_No_PLRIRR,sg_rep_grad_or_sanctions_No_PLRIRR_v, sg_rep_grad_or_sanctions_No_Q, sg_rep_grad_or_sanctions_No_Q_p, 
                               sg_rep_grad_or_sanctions_No_I2,sg_rep_grad_or_sanctions_grad_k, sg_rep_grad_or_sanctions_grad_PLRIRR,sg_rep_grad_or_sanctions_grad_PLRIRR_v, sg_rep_grad_or_sanctions_grad_Q, 
                               sg_rep_grad_or_sanctions_grad_Q_p, sg_rep_grad_or_sanctions_grad_I2, mod_rep_grad_or_sanctions_LRIRR, mod_rep_grad_or_sanctions_se,
                               mod_rep_grad_or_sanctions_res_df, mod_rep_grad_or_sanctions_Qres)

names(mod_rep_grad_or_sanctions_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                      "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                      "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_grad_or_sanctions_list)
########################################
### 8.4.19.2) Reported sanctions or graduated sanctions - Sanctions
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_grad_or_sanctions_sanct_names<- c('reported_sanctions_or_graduated_sanctions','programmatic', 'No', 'sanctions')
### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_grad_or_sanctions_No_es <-mods[mods$reported_sanctions_or_graduated_sanctions== 'No', "LRIRR"]
sg_rep_grad_or_sanctions_No_v <-mods[mods$reported_sanctions_or_graduated_sanctions== 'No', "RIRR_var"]
sg_rep_grad_or_sanctions_No_k <-length(sg_rep_grad_or_sanctions_No_es)
sg_rep_grad_or_sanctions_No_PLRIRR <-(sum(sg_rep_grad_or_sanctions_No_es/sg_rep_grad_or_sanctions_No_v))/(sum(1/sg_rep_grad_or_sanctions_No_v))
sg_rep_grad_or_sanctions_No_PLRIRR_v <- 1/sum(1/sg_rep_grad_or_sanctions_No_v)
sg_rep_grad_or_sanctions_No_Q <-sum((sg_rep_grad_or_sanctions_No_es-sg_rep_grad_or_sanctions_No_PLRIRR)^2/sg_rep_grad_or_sanctions_No_v)
sg_rep_grad_or_sanctions_No_Q_p <-pchisq(sg_rep_grad_or_sanctions_No_Q, df=(length(sg_rep_grad_or_sanctions_No_es)-1), lower.tail = FALSE)
sg_rep_grad_or_sanctions_No_I2 <-((sg_rep_grad_or_sanctions_No_Q-(length(sg_rep_grad_or_sanctions_No_es)-1))/sg_rep_grad_or_sanctions_No_Q)*100 

sg_rep_grad_or_sanctions_yes_es <-mods[mods$reported_sanctions_or_graduated_sanctions== 'Yes', "LRIRR"]
sg_rep_grad_or_sanctions_yes_v <-mods[mods$reported_sanctions_or_graduated_sanctions== 'Yes', "RIRR_var"]
sg_rep_grad_or_sanctions_yes_k <-length(sg_rep_grad_or_sanctions_yes_es)
sg_rep_grad_or_sanctions_yes_PLRIRR <-(sum(sg_rep_grad_or_sanctions_yes_es/sg_rep_grad_or_sanctions_yes_v))/(sum(1/sg_rep_grad_or_sanctions_yes_v))
sg_rep_grad_or_sanctions_yes_PLRIRR_v <- 1/sum(1/sg_rep_grad_or_sanctions_yes_v)
sg_rep_grad_or_sanctions_yes_Q <-sum((sg_rep_grad_or_sanctions_yes_es-sg_rep_grad_or_sanctions_yes_PLRIRR)^2/sg_rep_grad_or_sanctions_yes_v)
sg_rep_grad_or_sanctions_yes_Q_p <-pchisq(sg_rep_grad_or_sanctions_yes_Q, df=(length(sg_rep_grad_or_sanctions_yes_es)-1), lower.tail = FALSE)
sg_rep_grad_or_sanctions_yes_I2 <-((sg_rep_grad_or_sanctions_yes_Q-(length(sg_rep_grad_or_sanctions_yes_es)-1))/sg_rep_grad_or_sanctions_yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.grad.or.sanct<-lm(LRIRR ~ 1 + mods$reported_sanctions_or_graduated_sanctions,weights = 1/v)

###Extract moderator Beta
mod_rep_grad_or_sanctions_sanct_LRIRR<-summary(model.lm.mod.rep.grad.or.sanct)$coefficients[3, 1]

###Extract the Beta's standard error
mod_rep_grad_or_sanctions_sanct_se<-(summary(model.lm.mod.rep.grad.or.sanct)$coefficients[3, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_grad_or_sanctions_sanct_res_df<- model.lm.mod.rep.grad.or.sanct$df.residual

###Calculate Qres for the model
mod_rep_grad_or_sanctions_Y <- as.matrix(mods$LRIRR)
mod_rep_grad_or_sanctions_W <- diag(1/mods$RIRR_var)
mod_rep_grad_or_sanctions_X <- model.matrix(model.lm.mod.rep.grad.or.sanct)
mod_rep_grad_or_sanctions_P <- mod_rep_grad_or_sanctions_W - mod_rep_grad_or_sanctions_W %*% mod_rep_grad_or_sanctions_X %*% solve(
  t(mod_rep_grad_or_sanctions_X) %*% mod_rep_grad_or_sanctions_W %*% mod_rep_grad_or_sanctions_X) %*% t(
    mod_rep_grad_or_sanctions_X) %*% mod_rep_grad_or_sanctions_W 

mod_rep_grad_or_sanctions_Qres <- max(0, c(crossprod(mod_rep_grad_or_sanctions_Y, mod_rep_grad_or_sanctions_P) %*% mod_rep_grad_or_sanctions_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_grad_or_sanctions_sanct_list<-c(mod_rep_grad_or_sanctions_sanct_names,sg_rep_grad_or_sanctions_No_k, sg_rep_grad_or_sanctions_No_PLRIRR,sg_rep_grad_or_sanctions_No_PLRIRR_v, sg_rep_grad_or_sanctions_No_Q, sg_rep_grad_or_sanctions_No_Q_p, 
                                  sg_rep_grad_or_sanctions_No_I2,sg_rep_grad_or_sanctions_yes_k, sg_rep_grad_or_sanctions_yes_PLRIRR,sg_rep_grad_or_sanctions_yes_PLRIRR_v, sg_rep_grad_or_sanctions_yes_Q, 
                                  sg_rep_grad_or_sanctions_yes_Q_p, sg_rep_grad_or_sanctions_yes_I2, mod_rep_grad_or_sanctions_sanct_LRIRR, mod_rep_grad_or_sanctions_sanct_se,
                                  mod_rep_grad_or_sanctions_res_df, mod_rep_grad_or_sanctions_Qres)

names(mod_rep_grad_or_sanctions_sanct_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                         "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                         "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_grad_or_sanctions_sanct_list)
########################################
### 8.4.20) Reported using rewards
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_rewards_names<- c('reported_using_rewards','programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_rewards_No_es <-mods[mods$reported_using_rewards== 'No', "LRIRR"]
sg_rep_rewards_No_v <-mods[mods$reported_using_rewards== 'No', "RIRR_var"]
sg_rep_rewards_No_k <-length(sg_rep_rewards_No_es)
sg_rep_rewards_No_PLRIRR <-(sum(sg_rep_rewards_No_es/sg_rep_rewards_No_v))/(sum(1/sg_rep_rewards_No_v))
sg_rep_rewards_No_PLRIRR_v <- 1/sum(1/sg_rep_rewards_No_v)
sg_rep_rewards_No_Q <-sum((sg_rep_rewards_No_es-sg_rep_rewards_No_PLRIRR)^2/sg_rep_rewards_No_v)
sg_rep_rewards_No_Q_p <-pchisq(sg_rep_rewards_No_Q, df=(length(sg_rep_rewards_No_es)-1), lower.tail = FALSE)
sg_rep_rewards_No_I2 <-((sg_rep_rewards_No_Q-(length(sg_rep_rewards_No_es)-1))/sg_rep_rewards_No_Q)*100 

sg_rep_rewards_Yes_es <-mods[mods$reported_using_rewards== 'Yes', "LRIRR"]
sg_rep_rewards_Yes_v <-mods[mods$reported_using_rewards== 'Yes', "RIRR_var"]
sg_rep_rewards_Yes_k <-length(sg_rep_rewards_Yes_es)
sg_rep_rewards_Yes_PLRIRR <-(sum(sg_rep_rewards_Yes_es/sg_rep_rewards_Yes_v))/(sum(1/sg_rep_rewards_Yes_v))
sg_rep_rewards_Yes_PLRIRR_v <- 1/sum(1/sg_rep_rewards_Yes_v)
sg_rep_rewards_Yes_Q <-sum((sg_rep_rewards_Yes_es-sg_rep_rewards_Yes_PLRIRR)^2/sg_rep_rewards_Yes_v)
sg_rep_rewards_Yes_Q_p <-pchisq(sg_rep_rewards_Yes_Q, df=(length(sg_rep_rewards_Yes_es)-1), lower.tail = FALSE)
sg_rep_rewards_Yes_I2 <-((sg_rep_rewards_Yes_Q-(length(sg_rep_rewards_Yes_es)-1))/sg_rep_rewards_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.rewards<-lm(LRIRR ~ 1 + mods$reported_using_rewards,weights = 1/v)

###Extract moderator Beta
mod_rep_rewards_LRIRR<-summary(model.lm.mod.rep.rewards)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_rewards_se<-(summary(model.lm.mod.rep.rewards)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_rewards_res_df<- model.lm.mod.rep.rewards$df.residual

###Calculate Qres for the model
mod_rep_rewards_Y <- as.matrix(mods$LRIRR)
mod_rep_rewards_W <- diag(1/mods$RIRR_var)
mod_rep_rewards_X <- model.matrix(model.lm.mod.rep.rewards)
mod_rep_rewards_P <- mod_rep_rewards_W - mod_rep_rewards_W %*% mod_rep_rewards_X %*% solve(
  t(mod_rep_rewards_X) %*% mod_rep_rewards_W %*% mod_rep_rewards_X) %*% t(
    mod_rep_rewards_X) %*% mod_rep_rewards_W 

mod_rep_rewards_Qres <- max(0, c(crossprod(mod_rep_rewards_Y, mod_rep_rewards_P) %*% mod_rep_rewards_Y))
###Bind the list of values to the dataframe 'moderator_output'
mod_rep_rewards_list<-c(mod_rep_rewards_names,sg_rep_rewards_No_k, sg_rep_rewards_No_PLRIRR,sg_rep_rewards_No_PLRIRR_v, sg_rep_rewards_No_Q, sg_rep_rewards_No_Q_p, 
                             sg_rep_rewards_No_I2,sg_rep_rewards_Yes_k, sg_rep_rewards_Yes_PLRIRR,sg_rep_rewards_Yes_PLRIRR_v, sg_rep_rewards_Yes_Q, 
                             sg_rep_rewards_Yes_Q_p, sg_rep_rewards_Yes_I2, mod_rep_rewards_LRIRR, mod_rep_rewards_se,
                             mod_rep_rewards_res_df, mod_rep_rewards_Qres)

names(mod_rep_rewards_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                    "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                    "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_rewards_list)
########################################
### 8.4.21) Reported using sanctions and rewards
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_sanct_reward_names<- c('reported_using_sanctions_and_rewards','programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_sanct_reward_No_es <-mods[mods$reported_sanctions_and_rewards== 'No', "LRIRR"]
sg_rep_sanct_reward_No_v <-mods[mods$reported_sanctions_and_rewards== 'No', "RIRR_var"]
sg_rep_sanct_reward_No_k <-length(sg_rep_sanct_reward_No_es)
sg_rep_sanct_reward_No_PLRIRR <-(sum(sg_rep_sanct_reward_No_es/sg_rep_sanct_reward_No_v))/(sum(1/sg_rep_sanct_reward_No_v))
sg_rep_sanct_reward_No_PLRIRR_v <- 1/sum(1/sg_rep_sanct_reward_No_v)
sg_rep_sanct_reward_No_Q <-sum((sg_rep_sanct_reward_No_es-sg_rep_sanct_reward_No_PLRIRR)^2/sg_rep_sanct_reward_No_v)
sg_rep_sanct_reward_No_Q_p <-pchisq(sg_rep_sanct_reward_No_Q, df=(length(sg_rep_sanct_reward_No_es)-1), lower.tail = FALSE)
sg_rep_sanct_reward_No_I2 <-((sg_rep_sanct_reward_No_Q-(length(sg_rep_sanct_reward_No_es)-1))/sg_rep_sanct_reward_No_Q)*100 

sg_rep_sanct_reward_Yes_es <-mods[mods$reported_sanctions_and_rewards== 'Yes', "LRIRR"]
sg_rep_sanct_reward_Yes_v <-mods[mods$reported_sanctions_and_rewards== 'Yes', "RIRR_var"]
sg_rep_sanct_reward_Yes_k <-length(sg_rep_sanct_reward_Yes_es)
sg_rep_sanct_reward_Yes_PLRIRR <-(sum(sg_rep_sanct_reward_Yes_es/sg_rep_sanct_reward_Yes_v))/(sum(1/sg_rep_sanct_reward_Yes_v))
sg_rep_sanct_reward_Yes_PLRIRR_v <- 1/sum(1/sg_rep_sanct_reward_Yes_v)
sg_rep_sanct_reward_Yes_Q <-sum((sg_rep_sanct_reward_Yes_es-sg_rep_sanct_reward_Yes_PLRIRR)^2/sg_rep_sanct_reward_Yes_v)
sg_rep_sanct_reward_Yes_Q_p <-pchisq(sg_rep_sanct_reward_Yes_Q, df=(length(sg_rep_sanct_reward_Yes_es)-1), lower.tail = FALSE)
sg_rep_sanct_reward_Yes_I2 <-((sg_rep_sanct_reward_Yes_Q-(length(sg_rep_sanct_reward_Yes_es)-1))/sg_rep_sanct_reward_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.sanct.reward<-lm(LRIRR ~ 1 + mods$reported_sanctions_and_rewards,weights = 1/v)

###Extract moderator Beta
mod_rep_sanct_reward_LRIRR<-summary(model.lm.mod.rep.sanct.reward)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_sanct_reward_se<-(summary(model.lm.mod.rep.sanct.reward)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_sanct_reward_res_df<- model.lm.mod.rep.sanct.reward$df.residual

###Calculate Qres for the model
mod_rep_sanct_reward_Y <- as.matrix(mods$LRIRR)
mod_rep_sanct_reward_W <- diag(1/mods$RIRR_var)
mod_rep_sanct_reward_X <- model.matrix(model.lm.mod.rep.sanct.reward)
mod_rep_sanct_reward_P <- mod_rep_sanct_reward_W - mod_rep_sanct_reward_W %*% mod_rep_sanct_reward_X %*% solve(
  t(mod_rep_sanct_reward_X) %*% mod_rep_sanct_reward_W %*% mod_rep_sanct_reward_X) %*% t(
    mod_rep_sanct_reward_X) %*% mod_rep_sanct_reward_W 

mod_rep_sanct_reward_Qres <- max(0, c(crossprod(mod_rep_sanct_reward_Y, mod_rep_sanct_reward_P) %*% mod_rep_sanct_reward_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_sanct_reward_list<-c(mod_rep_sanct_reward_names,sg_rep_sanct_reward_No_k, sg_rep_sanct_reward_No_PLRIRR,sg_rep_sanct_reward_No_PLRIRR_v, sg_rep_sanct_reward_No_Q, sg_rep_sanct_reward_No_Q_p, 
                             sg_rep_sanct_reward_No_I2,sg_rep_sanct_reward_Yes_k, sg_rep_sanct_reward_Yes_PLRIRR,sg_rep_sanct_reward_Yes_PLRIRR_v, sg_rep_sanct_reward_Yes_Q, 
                             sg_rep_sanct_reward_Yes_Q_p, sg_rep_sanct_reward_Yes_I2, mod_rep_sanct_reward_LRIRR, mod_rep_sanct_reward_se,
                             mod_rep_sanct_reward_res_df, mod_rep_sanct_reward_Qres)

names(mod_rep_sanct_reward_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                    "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                    "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_sanct_reward_list)

########################################
### 8.4.22.1) Frequency of review hearings -  Weekly
########################################
### Create a a new variable, re-code '3 hearings per month' as NA because it has insufficient cells for moderation
mods$Frequency_review_hearing_recode <- na_if(mods$Frequency_review_hearing, '3 hearings per month')
if (class(mods$Frequency_review_hearing_recode ) == 'factor') {relevel(mods$Frequency_review_hearing_recode, ref = "fortnightly")} # re-level if need be

### Create an object for the row name, level_ref, and level_alternative names
mod_freq_rh_weekly_names<- c('frequency_of_review_hearings', 'programmatic', 'Fortnightly', 'Weekly')

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_freq_rh_fortnightly_es <-na.omit(mods[mods$Frequency_review_hearing_recode == 'fortnightly', "LRIRR"])
sg_freq_rh_fortnightly_v <-na.omit(mods[mods$Frequency_review_hearing_recode == 'fortnightly', "RIRR_var"])
sg_freq_rh_fortnightly_k <-length(sg_freq_rh_fortnightly_es)
sg_freq_rh_fortnightly_PLRIRR <-(sum(sg_freq_rh_fortnightly_es/sg_freq_rh_fortnightly_v))/(sum(1/sg_freq_rh_fortnightly_v))
sg_freq_rh_fortnightly_PLRIRR_v <- 1/sum(1/sg_freq_rh_fortnightly_v)
sg_freq_rh_fortnightly_Q <-sum((sg_freq_rh_fortnightly_es-sg_freq_rh_fortnightly_PLRIRR)^2/sg_freq_rh_fortnightly_v)
sg_freq_rh_fortnightly_Q_p <-pchisq(sg_freq_rh_fortnightly_Q, df=(length(sg_freq_rh_fortnightly_es)-1), lower.tail = FALSE)
sg_freq_rh_fortnightly_I2 <-((sg_freq_rh_fortnightly_Q-(length(sg_freq_rh_fortnightly_es)-1))/sg_freq_rh_fortnightly_Q)*100

sg_freq_rh_weekly_es <-na.omit(mods[mods$Frequency_review_hearing_recode == 'Weekly', "LRIRR"])
sg_freq_rh_weekly_v <-na.omit(mods[mods$Frequency_review_hearing_recode == 'Weekly', "RIRR_var"])
sg_freq_rh_weekly_k <-length(sg_freq_rh_weekly_es)
sg_freq_rh_weekly_PLRIRR <-(sum(sg_freq_rh_weekly_es/sg_freq_rh_weekly_v))/(sum(1/sg_freq_rh_weekly_v))
sg_freq_rh_weekly_PLRIRR_v <- 1/sum(1/sg_freq_rh_weekly_v)
sg_freq_rh_weekly_Q <-sum((sg_freq_rh_weekly_es-sg_freq_rh_weekly_PLRIRR)^2/sg_freq_rh_weekly_v)
sg_freq_rh_weekly_Q_p <-pchisq(sg_freq_rh_weekly_Q, df=(length(sg_freq_rh_weekly_es)-1), lower.tail = FALSE)
sg_freq_rh_weekly_I2 <-((sg_freq_rh_weekly_Q-(length(sg_freq_rh_weekly_es)-1))/sg_freq_rh_weekly_Q)*100

###Insert the moderator into the model
model.lm.mod.freq.rh<-lm(LRIRR ~ 1 + mods$Frequency_review_hearing_recode,weights = 1/v)

###Extract moderator Beta
mod_freq_rh_weekly_LRIRR<-summary(model.lm.mod.freq.rh)$coefficients[3, 1]

###Extract the Beta's standard error of the model
mod_freq_rh_weekly_se<-summary(model.lm.mod.freq.rh)$coefficients[3, 2]

###Create an object for the residual degrees of freedom. 
mod_freq_rh_weekly_res_df<- model.lm.mod.freq.rh$df.residual

###Calculate Qres for the model
mod_freq_rh_es <-mods[mods$Frequency_review_hearing_recode %in% na.omit(mods$Frequency_review_hearing_recode), "LRIRR"]
mod_freq_rh_v <-mods[mods$Frequency_review_hearing_recode %in% na.omit(mods$Frequency_review_hearing_recode), "RIRR_var"]

mod_freq_rh_Y <- as.matrix(mod_freq_rh_es)
mod_freq_rh_W <- diag(1/mod_freq_rh_v)
mod_freq_rh_X <- model.matrix(model.lm.mod.freq.rh)
mod_freq_rh_P <- mod_freq_rh_W - mod_freq_rh_W %*% mod_freq_rh_X %*% solve(
  t(mod_freq_rh_X) %*% mod_freq_rh_W %*% mod_freq_rh_X) %*% t(
    mod_freq_rh_X) %*% mod_freq_rh_W 

mod_freq_rh_Qres <- max(0, c(crossprod(mod_freq_rh_Y, mod_freq_rh_P) %*% mod_freq_rh_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_freq_rh_weekly_list<-c(mod_freq_rh_weekly_names,sg_freq_rh_fortnightly_k, sg_freq_rh_fortnightly_PLRIRR,sg_freq_rh_fortnightly_PLRIRR_v, sg_freq_rh_fortnightly_Q, sg_freq_rh_fortnightly_Q_p, 
                              sg_freq_rh_fortnightly_I2,sg_freq_rh_weekly_k, sg_freq_rh_weekly_PLRIRR,sg_freq_rh_weekly_PLRIRR_v, sg_freq_rh_weekly_Q, sg_freq_rh_weekly_Q_p, 
                              sg_freq_rh_weekly_I2, mod_freq_rh_weekly_LRIRR, mod_freq_rh_weekly_se, mod_freq_rh_weekly_res_df, 
                              mod_freq_rh_Qres)

names(mod_freq_rh_weekly_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                     "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                     "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_freq_rh_weekly_list)

########################################
### 8.4.22.2) Frequency of review hearings -  Matched to level of risk and need
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_freq_rh_matched_to_rn_names<- c('frequency_of_review_hearings', 'programmatic', 'Fortnightly', 'Matched to risk and need')

### Calculate the PLRIRR, Q, P value for Q, and I^2 for each subgroup
sg_freq_rh_fortnightly_es <-na.omit(mods[mods$Frequency_review_hearing_recode == 'fortnightly', "LRIRR"])
sg_freq_rh_fortnightly_v <-na.omit(mods[mods$Frequency_review_hearing_recode == 'fortnightly', "RIRR_var"])
sg_freq_rh_fortnightly_k <-length(sg_freq_rh_fortnightly_es)
sg_freq_rh_fortnightly_PLRIRR <-(sum(sg_freq_rh_fortnightly_es/sg_freq_rh_fortnightly_v))/(sum(1/sg_freq_rh_fortnightly_v))
sg_freq_rh_fortnightly_PLRIRR_v <- 1/sum(1/sg_freq_rh_fortnightly_v)
sg_freq_rh_fortnightly_Q <-sum((sg_freq_rh_fortnightly_es-sg_freq_rh_fortnightly_PLRIRR)^2/sg_freq_rh_fortnightly_v)
sg_freq_rh_fortnightly_Q_p <-pchisq(sg_freq_rh_fortnightly_Q, df=(length(sg_freq_rh_fortnightly_es)-1), lower.tail = FALSE)
sg_freq_rh_fortnightly_I2 <-((sg_freq_rh_fortnightly_Q-(length(sg_freq_rh_fortnightly_es)-1))/sg_freq_rh_fortnightly_Q)*100

sg_freq_rh_matched_to_rn_es <-na.omit(mods[mods$Frequency_review_hearing_recode == 'matched_to_rn', "LRIRR"])
sg_freq_rh_matched_to_rn_v <-na.omit(mods[mods$Frequency_review_hearing_recode == 'matched_to_rn', "RIRR_var"])
sg_freq_rh_matched_to_rn_k <-length(sg_freq_rh_matched_to_rn_es)
sg_freq_rh_matched_to_rn_PLRIRR <-(sum(sg_freq_rh_matched_to_rn_es/sg_freq_rh_matched_to_rn_v))/(sum(1/sg_freq_rh_matched_to_rn_v))
sg_freq_rh_matched_to_rn_PLRIRR_v <- 1/sum(1/sg_freq_rh_matched_to_rn_v)
sg_freq_rh_matched_to_rn_Q <-sum((sg_freq_rh_matched_to_rn_es-sg_freq_rh_matched_to_rn_PLRIRR)^2/sg_freq_rh_matched_to_rn_v)
sg_freq_rh_matched_to_rn_Q_p <-pchisq(sg_freq_rh_matched_to_rn_Q, df=(length(sg_freq_rh_matched_to_rn_es)-1), lower.tail = FALSE)
sg_freq_rh_matched_to_rn_I2 <-((sg_freq_rh_matched_to_rn_Q-(length(sg_freq_rh_matched_to_rn_es)-1))/sg_freq_rh_matched_to_rn_Q)*100

###Insert the moderator into the model
model.lm.mod.freq.rh<-lm(LRIRR ~ 1 + mods$Frequency_review_hearing_recode,weights = 1/v)

###Extract moderator Beta
mod_freq_rh_matched_to_rn_LRIRR<-summary(model.lm.mod.freq.rh)$coefficients[2, 1]

###Extract the Beta's standard error of the model
mod_freq_rh_matched_to_rn_se<-summary(model.lm.mod.freq.rh)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_freq_rh_matched_to_rn_res_df<- model.lm.mod.freq.rh$df.residual

###Calculate Qres for the model
mod_freq_rh_es <-mods[mods$Frequency_review_hearing_recode %in% na.omit(mods$Frequency_review_hearing_recode), "LRIRR"]
mod_freq_rh_v <-mods[mods$Frequency_review_hearing_recode %in% na.omit(mods$Frequency_review_hearing_recode), "RIRR_var"]

mod_freq_rh_Y <- as.matrix(mod_freq_rh_es)
mod_freq_rh_W <- diag(1/mod_freq_rh_v)
mod_freq_rh_X <- model.matrix(model.lm.mod.freq.rh)
mod_freq_rh_P <- mod_freq_rh_W - mod_freq_rh_W %*% mod_freq_rh_X %*% solve(
  t(mod_freq_rh_X) %*% mod_freq_rh_W %*% mod_freq_rh_X) %*% t(
    mod_freq_rh_X) %*% mod_freq_rh_W 

mod_freq_rh_Qres <- max(0, c(crossprod(mod_freq_rh_Y, mod_freq_rh_P) %*% mod_freq_rh_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_freq_rh_matched_to_rn_list<-c(mod_freq_rh_matched_to_rn_names,sg_freq_rh_fortnightly_k, sg_freq_rh_fortnightly_PLRIRR,sg_freq_rh_fortnightly_PLRIRR_v, sg_freq_rh_fortnightly_Q, sg_freq_rh_fortnightly_Q_p, 
                                 sg_freq_rh_fortnightly_I2,sg_freq_rh_matched_to_rn_k, sg_freq_rh_matched_to_rn_PLRIRR,sg_freq_rh_matched_to_rn_PLRIRR_v, sg_freq_rh_matched_to_rn_Q, sg_freq_rh_matched_to_rn_Q_p, 
                                 sg_freq_rh_matched_to_rn_I2, mod_freq_rh_matched_to_rn_LRIRR, mod_freq_rh_matched_to_rn_se, mod_freq_rh_matched_to_rn_res_df, 
                                 mod_freq_rh_Qres)

names(mod_freq_rh_matched_to_rn_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                        "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                        "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_freq_rh_matched_to_rn_list)

########################################
### 8.4.23) Frequent review hearings collapsed
########################################
###Re-level the variable if need be 
if (class(mods$frequent_rh_collapsed ) == 'factor') {relevel(mods$frequent_rh_collapsed , ref = "No")} # re-level if need be

### Create an object for the row name, level_ref, and level_alternative names
mod_rep_freq_rh_names <-c('reported_frequent_review_hearings','programmatic', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_freq_rh_No_es <-mods[mods$frequent_rh_collapsed== 'No', "LRIRR"]
sg_rep_freq_rh_No_v <-mods[mods$frequent_rh_collapsed== 'No', "RIRR_var"]
sg_rep_freq_rh_No_k <-length(sg_rep_freq_rh_No_es)
sg_rep_freq_rh_No_PLRIRR <-(sum(sg_rep_freq_rh_No_es/sg_rep_freq_rh_No_v))/(sum(1/sg_rep_freq_rh_No_v))
sg_rep_freq_rh_No_PLRIRR_v <- 1/sum(1/sg_rep_freq_rh_No_v)
sg_rep_freq_rh_No_Q <-sum((sg_rep_freq_rh_No_es-sg_rep_freq_rh_No_PLRIRR)^2/sg_rep_freq_rh_No_v)
sg_rep_freq_rh_No_Q_p <-pchisq(sg_rep_freq_rh_No_Q, df=(length(sg_rep_freq_rh_No_es)-1), lower.tail = FALSE)
sg_rep_freq_rh_No_I2 <-((sg_rep_freq_rh_No_Q-(length(sg_rep_freq_rh_No_es)-1))/sg_rep_freq_rh_No_Q)*100 

sg_rep_freq_rh_Yes_es <-mods[mods$frequent_rh_collapsed== 'frequent', "LRIRR"]
sg_rep_freq_rh_Yes_v <-mods[mods$frequent_rh_collapsed== 'frequent', "RIRR_var"]
sg_rep_freq_rh_Yes_k <-length(sg_rep_freq_rh_Yes_es)
sg_rep_freq_rh_Yes_PLRIRR <-(sum(sg_rep_freq_rh_Yes_es/sg_rep_freq_rh_Yes_v))/(sum(1/sg_rep_freq_rh_Yes_v))
sg_rep_freq_rh_Yes_PLRIRR_v <- 1/sum(1/sg_rep_freq_rh_Yes_v)
sg_rep_freq_rh_Yes_Q <-sum((sg_rep_freq_rh_Yes_es-sg_rep_freq_rh_Yes_PLRIRR)^2/sg_rep_freq_rh_Yes_v)
sg_rep_freq_rh_Yes_Q_p <-pchisq(sg_rep_freq_rh_Yes_Q, df=(length(sg_rep_freq_rh_Yes_es)-1), lower.tail = FALSE)
sg_rep_freq_rh_Yes_I2 <-((sg_rep_freq_rh_Yes_Q-(length(sg_rep_freq_rh_Yes_es)-1))/sg_rep_freq_rh_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.freq.rh<-lm(LRIRR ~ 1 + mods$frequent_rh_collapsed,weights = 1/v)

###Extract moderator Beta
mod_rep_freq_rh_LRIRR<-summary(model.lm.mod.rep.freq.rh)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_freq_rh_se<-(summary(model.lm.mod.rep.freq.rh)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_freq_rh_res_df<- model.lm.mod.rep.freq.rh$df.residual

###Calculate Qres for the model
mod_rep_freq_rh_Y <- as.matrix(mods$LRIRR)
mod_rep_freq_rh_W <- diag(1/mods$RIRR_var)
mod_rep_freq_rh_X <- model.matrix(model.lm.mod.rep.freq.rh)
mod_rep_freq_rh_P <- mod_rep_freq_rh_W - mod_rep_freq_rh_W %*% mod_rep_freq_rh_X %*% solve(
  t(mod_rep_freq_rh_X) %*% mod_rep_freq_rh_W %*% mod_rep_freq_rh_X) %*% t(
    mod_rep_freq_rh_X) %*% mod_rep_freq_rh_W 

mod_rep_freq_rh_Qres <- max(0, c(crossprod(mod_rep_freq_rh_Y, mod_rep_freq_rh_P) %*% mod_rep_freq_rh_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_freq_rh_list<-c(mod_rep_freq_rh_names,sg_rep_freq_rh_No_k, sg_rep_freq_rh_No_PLRIRR,sg_rep_freq_rh_No_PLRIRR_v, sg_rep_freq_rh_No_Q, sg_rep_freq_rh_No_Q_p, 
                        sg_rep_freq_rh_No_I2,sg_rep_freq_rh_Yes_k, sg_rep_freq_rh_Yes_PLRIRR,sg_rep_freq_rh_Yes_PLRIRR_v, sg_rep_freq_rh_Yes_Q, 
                        sg_rep_freq_rh_Yes_Q_p, sg_rep_freq_rh_Yes_I2, mod_rep_freq_rh_LRIRR, mod_rep_freq_rh_se,
                        mod_rep_freq_rh_res_df, mod_rep_freq_rh_Qres)

names(mod_rep_freq_rh_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                               "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                               "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_freq_rh_list)
########################################
### 8.5) Treatment
########################################
########################################
### 8.5.1) Reported measuring risk - Treatment
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_meas_risk_tr_names <-c('reported_measuring_risk_treat','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_meas_risk_tr_No_es <-mods[mods$reported_measuring_risk_treat== 'No', "LRIRR"]
sg_rep_meas_risk_tr_No_v <-mods[mods$reported_measuring_risk_treat== 'No', "RIRR_var"]
sg_rep_meas_risk_tr_No_k <-length(sg_rep_meas_risk_tr_No_es)
sg_rep_meas_risk_tr_No_PLRIRR <-(sum(sg_rep_meas_risk_tr_No_es/sg_rep_meas_risk_tr_No_v))/(sum(1/sg_rep_meas_risk_tr_No_v))
sg_rep_meas_risk_tr_No_PLRIRR_v <- 1/sum(1/sg_rep_meas_risk_tr_No_v)
sg_rep_meas_risk_tr_No_Q <-sum((sg_rep_meas_risk_tr_No_es-sg_rep_meas_risk_tr_No_PLRIRR)^2/sg_rep_meas_risk_tr_No_v)
sg_rep_meas_risk_tr_No_Q_p <-pchisq(sg_rep_meas_risk_tr_No_Q, df=(length(sg_rep_meas_risk_tr_No_es)-1), lower.tail = FALSE)
sg_rep_meas_risk_tr_No_I2 <-((sg_rep_meas_risk_tr_No_Q-(length(sg_rep_meas_risk_tr_No_es)-1))/sg_rep_meas_risk_tr_No_Q)*100 

sg_rep_meas_risk_tr_Yes_es <-mods[mods$reported_measuring_risk_treat== 'Yes', "LRIRR"]
sg_rep_meas_risk_tr_Yes_v <-mods[mods$reported_measuring_risk_treat== 'Yes', "RIRR_var"]
sg_rep_meas_risk_tr_Yes_k <-length(sg_rep_meas_risk_tr_Yes_es)
sg_rep_meas_risk_tr_Yes_PLRIRR <-(sum(sg_rep_meas_risk_tr_Yes_es/sg_rep_meas_risk_tr_Yes_v))/(sum(1/sg_rep_meas_risk_tr_Yes_v))
sg_rep_meas_risk_tr_Yes_PLRIRR_v <- 1/sum(1/sg_rep_meas_risk_tr_Yes_v)
sg_rep_meas_risk_tr_Yes_Q <-sum((sg_rep_meas_risk_tr_Yes_es-sg_rep_meas_risk_tr_Yes_PLRIRR)^2/sg_rep_meas_risk_tr_Yes_v)
sg_rep_meas_risk_tr_Yes_Q_p <-pchisq(sg_rep_meas_risk_tr_Yes_Q, df=(length(sg_rep_meas_risk_tr_Yes_es)-1), lower.tail = FALSE)
sg_rep_meas_risk_tr_Yes_I2 <-((sg_rep_meas_risk_tr_Yes_Q-(length(sg_rep_meas_risk_tr_Yes_es)-1))/sg_rep_meas_risk_tr_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.meas.risk.tr<-lm(LRIRR ~ 1 + mods$reported_measuring_risk_treat,weights = 1/v)

###Extract moderator Beta
mod_rep_meas_risk_tr_LRIRR<-summary(model.lm.mod.rep.meas.risk.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_meas_risk_tr_se<-(summary(model.lm.mod.rep.meas.risk.tr)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_meas_risk_tr_res_df<- model.lm.mod.rep.meas.risk.tr$df.residual

###Calculate Qres for the model
mod_rep_meas_risk_tr_Y <- as.matrix(mods$LRIRR)
mod_rep_meas_risk_tr_W <- diag(1/mods$RIRR_var)
mod_rep_meas_risk_tr_X <- model.matrix(model.lm.mod.rep.meas.risk.tr)
mod_rep_meas_risk_tr_P <- mod_rep_meas_risk_tr_W - mod_rep_meas_risk_tr_W %*% mod_rep_meas_risk_tr_X %*% solve(
  t(mod_rep_meas_risk_tr_X) %*% mod_rep_meas_risk_tr_W %*% mod_rep_meas_risk_tr_X) %*% t(
    mod_rep_meas_risk_tr_X) %*% mod_rep_meas_risk_tr_W 

mod_rep_meas_risk_tr_Qres <- max(0, c(crossprod(mod_rep_meas_risk_tr_Y, mod_rep_meas_risk_tr_P) %*% mod_rep_meas_risk_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_meas_risk_tr_list<-c(mod_rep_meas_risk_tr_names,sg_rep_meas_risk_tr_No_k, sg_rep_meas_risk_tr_No_PLRIRR,sg_rep_meas_risk_tr_No_PLRIRR_v, sg_rep_meas_risk_tr_No_Q, sg_rep_meas_risk_tr_No_Q_p, 
                             sg_rep_meas_risk_tr_No_I2,sg_rep_meas_risk_tr_Yes_k, sg_rep_meas_risk_tr_Yes_PLRIRR,sg_rep_meas_risk_tr_Yes_PLRIRR_v, sg_rep_meas_risk_tr_Yes_Q, 
                             sg_rep_meas_risk_tr_Yes_Q_p, sg_rep_meas_risk_tr_Yes_I2, mod_rep_meas_risk_tr_LRIRR, mod_rep_meas_risk_tr_se,
                             mod_rep_meas_risk_tr_res_df, mod_rep_meas_risk_tr_Qres)

names(mod_rep_meas_risk_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                    "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                    "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_meas_risk_tr_list)
########################################
### 8.5.2) Reported measuring risk - Control
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_meas_risk_con_names <-c('reported_measuring_risk_cont','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_meas_risk_con_No_es <-mods[mods$reported_measuring_risk_cont== 'No', "LRIRR"]
sg_rep_meas_risk_con_No_v <-mods[mods$reported_measuring_risk_cont== 'No', "RIRR_var"]
sg_rep_meas_risk_con_No_k <-length(sg_rep_meas_risk_con_No_es)
sg_rep_meas_risk_con_No_PLRIRR <-(sum(sg_rep_meas_risk_con_No_es/sg_rep_meas_risk_con_No_v))/(sum(1/sg_rep_meas_risk_con_No_v))
sg_rep_meas_risk_con_No_PLRIRR_v <- 1/sum(1/sg_rep_meas_risk_con_No_v)
sg_rep_meas_risk_con_No_Q <-sum((sg_rep_meas_risk_con_No_es-sg_rep_meas_risk_con_No_PLRIRR)^2/sg_rep_meas_risk_con_No_v)
sg_rep_meas_risk_con_No_Q_p <-pchisq(sg_rep_meas_risk_con_No_Q, df=(length(sg_rep_meas_risk_con_No_es)-1), lower.tail = FALSE)
sg_rep_meas_risk_con_No_I2 <-((sg_rep_meas_risk_con_No_Q-(length(sg_rep_meas_risk_con_No_es)-1))/sg_rep_meas_risk_con_No_Q)*100 

sg_rep_meas_risk_con_Yes_es <-mods[mods$reported_measuring_risk_cont== 'Yes', "LRIRR"]
sg_rep_meas_risk_con_Yes_v <-mods[mods$reported_measuring_risk_cont== 'Yes', "RIRR_var"]
sg_rep_meas_risk_con_Yes_k <-length(sg_rep_meas_risk_con_Yes_es)
sg_rep_meas_risk_con_Yes_PLRIRR <-(sum(sg_rep_meas_risk_con_Yes_es/sg_rep_meas_risk_con_Yes_v))/(sum(1/sg_rep_meas_risk_con_Yes_v))
sg_rep_meas_risk_con_Yes_PLRIRR_v <- 1/sum(1/sg_rep_meas_risk_con_Yes_v)
sg_rep_meas_risk_con_Yes_Q <-sum((sg_rep_meas_risk_con_Yes_es-sg_rep_meas_risk_con_Yes_PLRIRR)^2/sg_rep_meas_risk_con_Yes_v)
sg_rep_meas_risk_con_Yes_Q_p <-pchisq(sg_rep_meas_risk_con_Yes_Q, df=(length(sg_rep_meas_risk_con_Yes_es)-1), lower.tail = FALSE)
sg_rep_meas_risk_con_Yes_I2 <-((sg_rep_meas_risk_con_Yes_Q-(length(sg_rep_meas_risk_con_Yes_es)-1))/sg_rep_meas_risk_con_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.meas.risk.con<-lm(LRIRR ~ 1 + mods$reported_measuring_risk_cont,weights = 1/v)

###Extract moderator Beta
mod_rep_meas_risk_con_LRIRR<-summary(model.lm.mod.rep.meas.risk.con)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_meas_risk_con_se<- summary(model.lm.mod.rep.meas.risk.con)$coefficients[2, 2] 

###Create an object for the residual degrees of freedom. 
mod_rep_meas_risk_con_res_df<- model.lm.mod.rep.meas.risk.con$df.residual

###Calculate Qres for the model
mod_rep_meas_risk_con_Y <- as.matrix(mods$LRIRR)
mod_rep_meas_risk_con_W <- diag(1/mods$RIRR_var)
mod_rep_meas_risk_con_X <- model.matrix(model.lm.mod.rep.meas.risk.con)
mod_rep_meas_risk_con_P <- mod_rep_meas_risk_con_W - mod_rep_meas_risk_con_W %*% mod_rep_meas_risk_con_X %*% solve(
  t(mod_rep_meas_risk_con_X) %*% mod_rep_meas_risk_con_W %*% mod_rep_meas_risk_con_X) %*% t(
    mod_rep_meas_risk_con_X) %*% mod_rep_meas_risk_con_W 

mod_rep_meas_risk_con_Qres <- max(0, c(crossprod(mod_rep_meas_risk_con_Y, mod_rep_meas_risk_con_P) %*% mod_rep_meas_risk_con_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_meas_risk_con_list<-c(mod_rep_meas_risk_con_names,sg_rep_meas_risk_con_No_k, sg_rep_meas_risk_con_No_PLRIRR,sg_rep_meas_risk_con_No_PLRIRR_v, sg_rep_meas_risk_con_No_Q, sg_rep_meas_risk_con_No_Q_p, 
                              sg_rep_meas_risk_con_No_I2,sg_rep_meas_risk_con_Yes_k, sg_rep_meas_risk_con_Yes_PLRIRR,sg_rep_meas_risk_con_Yes_PLRIRR_v, sg_rep_meas_risk_con_Yes_Q, 
                              sg_rep_meas_risk_con_Yes_Q_p, sg_rep_meas_risk_con_Yes_I2, mod_rep_meas_risk_con_LRIRR, mod_rep_meas_risk_con_se,
                              mod_rep_meas_risk_con_res_df, mod_rep_meas_risk_con_Qres)

names(mod_rep_meas_risk_con_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                     "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                     "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_meas_risk_con_list)
########################################
### 8.5.3) Reported individualising treatment - Treatment
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_indi_treat_tr_names <-c('reported_individualized_treat_treat','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_indi_treat_tr_No_es <-mods[mods$reported_individualized_treat_treat== 'No', "LRIRR"]
sg_rep_indi_treat_tr_No_v <-mods[mods$reported_individualized_treat_treat== 'No', "RIRR_var"]
sg_rep_indi_treat_tr_No_k <-length(sg_rep_indi_treat_tr_No_es)
sg_rep_indi_treat_tr_No_PLRIRR <-(sum(sg_rep_indi_treat_tr_No_es/sg_rep_indi_treat_tr_No_v))/(sum(1/sg_rep_indi_treat_tr_No_v))
sg_rep_indi_treat_tr_No_PLRIRR_v <- 1/sum(1/sg_rep_indi_treat_tr_No_v)
sg_rep_indi_treat_tr_No_Q <-sum((sg_rep_indi_treat_tr_No_es-sg_rep_indi_treat_tr_No_PLRIRR)^2/sg_rep_indi_treat_tr_No_v)
sg_rep_indi_treat_tr_No_Q_p <-pchisq(sg_rep_indi_treat_tr_No_Q, df=(length(sg_rep_indi_treat_tr_No_es)-1), lower.tail = FALSE)
sg_rep_indi_treat_tr_No_I2 <-((sg_rep_indi_treat_tr_No_Q-(length(sg_rep_indi_treat_tr_No_es)-1))/sg_rep_indi_treat_tr_No_Q)*100 

sg_rep_indi_treat_tr_Yes_es <-mods[mods$reported_individualized_treat_treat== 'Yes', "LRIRR"]
sg_rep_indi_treat_tr_Yes_v <-mods[mods$reported_individualized_treat_treat== 'Yes', "RIRR_var"]
sg_rep_indi_treat_tr_Yes_k <-length(sg_rep_indi_treat_tr_Yes_es)
sg_rep_indi_treat_tr_Yes_PLRIRR <-(sum(sg_rep_indi_treat_tr_Yes_es/sg_rep_indi_treat_tr_Yes_v))/(sum(1/sg_rep_indi_treat_tr_Yes_v))
sg_rep_indi_treat_tr_Yes_PLRIRR_v <- 1/sum(1/sg_rep_indi_treat_tr_Yes_v)
sg_rep_indi_treat_tr_Yes_Q <-sum((sg_rep_indi_treat_tr_Yes_es-sg_rep_indi_treat_tr_Yes_PLRIRR)^2/sg_rep_indi_treat_tr_Yes_v)
sg_rep_indi_treat_tr_Yes_Q_p <-pchisq(sg_rep_indi_treat_tr_Yes_Q, df=(length(sg_rep_indi_treat_tr_Yes_es)-1), lower.tail = FALSE)
sg_rep_indi_treat_tr_Yes_I2 <-((sg_rep_indi_treat_tr_Yes_Q-(length(sg_rep_indi_treat_tr_Yes_es)-1))/sg_rep_indi_treat_tr_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.indi.treat.tr<-lm(LRIRR ~ 1 + mods$reported_individualized_treat_treat,weights = 1/v)

###Extract moderator Beta
mod_rep_indi_treat_tr_LRIRR<-summary(model.lm.mod.rep.indi.treat.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_indi_treat_tr_se<-(summary(model.lm.mod.rep.indi.treat.tr)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_indi_treat_tr_res_df<- model.lm.mod.rep.indi.treat.tr$df.residual

###Calculate Qres for the model
mod_rep_indi_treat_tr_Y <- as.matrix(mods$LRIRR)
mod_rep_indi_treat_tr_W <- diag(1/mods$RIRR_var)
mod_rep_indi_treat_tr_X <- model.matrix(model.lm.mod.rep.indi.treat.tr)
mod_rep_indi_treat_tr_P <- mod_rep_indi_treat_tr_W - mod_rep_indi_treat_tr_W %*% mod_rep_indi_treat_tr_X %*% solve(
  t(mod_rep_indi_treat_tr_X) %*% mod_rep_indi_treat_tr_W %*% mod_rep_indi_treat_tr_X) %*% t(
    mod_rep_indi_treat_tr_X) %*% mod_rep_indi_treat_tr_W 

mod_rep_indi_treat_tr_Qres <- max(0, c(crossprod(mod_rep_indi_treat_tr_Y, mod_rep_indi_treat_tr_P) %*% mod_rep_indi_treat_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_indi_treat_tr_list<-c(mod_rep_indi_treat_tr_names,sg_rep_indi_treat_tr_No_k, sg_rep_indi_treat_tr_No_PLRIRR,sg_rep_indi_treat_tr_No_PLRIRR_v, sg_rep_indi_treat_tr_No_Q, sg_rep_indi_treat_tr_No_Q_p, 
                              sg_rep_indi_treat_tr_No_I2,sg_rep_indi_treat_tr_Yes_k, sg_rep_indi_treat_tr_Yes_PLRIRR,sg_rep_indi_treat_tr_Yes_PLRIRR_v, sg_rep_indi_treat_tr_Yes_Q, 
                              sg_rep_indi_treat_tr_Yes_Q_p, sg_rep_indi_treat_tr_Yes_I2, mod_rep_indi_treat_tr_LRIRR, mod_rep_indi_treat_tr_se,
                              mod_rep_indi_treat_tr_res_df, mod_rep_indi_treat_tr_Qres)

names(mod_rep_indi_treat_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                     "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                     "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_indi_treat_tr_list)
########################################
### 8.5.4) Reported individualising treatment - Control
########################################
### Note: there was insufficient variance in this moderator to perform an analysis 
########################################
### 8.5.5) Reported drug or alcohol treatment - Treatment
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_drug_alco_treat_tr_names <-c('reported_drug_or_alcohol_treatment_treat','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_drug_alco_treat_tr_No_es <-mods[mods$reported_drug_or_alcohol_treatment_treat== 'No', "LRIRR"]
sg_rep_drug_alco_treat_tr_No_v <-mods[mods$reported_drug_or_alcohol_treatment_treat== 'No', "RIRR_var"]
sg_rep_drug_alco_treat_tr_No_k <-length(sg_rep_drug_alco_treat_tr_No_es)
sg_rep_drug_alco_treat_tr_No_PLRIRR <-(sum(sg_rep_drug_alco_treat_tr_No_es/sg_rep_drug_alco_treat_tr_No_v))/(sum(1/sg_rep_drug_alco_treat_tr_No_v))
sg_rep_drug_alco_treat_tr_No_PLRIRR_v <- 1/sum(1/sg_rep_drug_alco_treat_tr_No_v)
sg_rep_drug_alco_treat_tr_No_Q <-sum((sg_rep_drug_alco_treat_tr_No_es-sg_rep_drug_alco_treat_tr_No_PLRIRR)^2/sg_rep_drug_alco_treat_tr_No_v)
sg_rep_drug_alco_treat_tr_No_Q_p <-pchisq(sg_rep_drug_alco_treat_tr_No_Q, df=(length(sg_rep_drug_alco_treat_tr_No_es)-1), lower.tail = FALSE)
sg_rep_drug_alco_treat_tr_No_I2 <-((sg_rep_drug_alco_treat_tr_No_Q-(length(sg_rep_drug_alco_treat_tr_No_es)-1))/sg_rep_drug_alco_treat_tr_No_Q)*100 

sg_rep_drug_alco_treat_tr_Yes_es <-mods[mods$reported_drug_or_alcohol_treatment_treat== 'Yes', "LRIRR"]
sg_rep_drug_alco_treat_tr_Yes_v <-mods[mods$reported_drug_or_alcohol_treatment_treat== 'Yes', "RIRR_var"]
sg_rep_drug_alco_treat_tr_Yes_k <-length(sg_rep_drug_alco_treat_tr_Yes_es)
sg_rep_drug_alco_treat_tr_Yes_PLRIRR <-(sum(sg_rep_drug_alco_treat_tr_Yes_es/sg_rep_drug_alco_treat_tr_Yes_v))/(sum(1/sg_rep_drug_alco_treat_tr_Yes_v))
sg_rep_drug_alco_treat_tr_Yes_PLRIRR_v <- 1/sum(1/sg_rep_drug_alco_treat_tr_Yes_v)
sg_rep_drug_alco_treat_tr_Yes_Q <-sum((sg_rep_drug_alco_treat_tr_Yes_es-sg_rep_drug_alco_treat_tr_Yes_PLRIRR)^2/sg_rep_drug_alco_treat_tr_Yes_v)
sg_rep_drug_alco_treat_tr_Yes_Q_p <-pchisq(sg_rep_drug_alco_treat_tr_Yes_Q, df=(length(sg_rep_drug_alco_treat_tr_Yes_es)-1), lower.tail = FALSE)
sg_rep_drug_alco_treat_tr_Yes_I2 <-((sg_rep_drug_alco_treat_tr_Yes_Q-(length(sg_rep_drug_alco_treat_tr_Yes_es)-1))/sg_rep_drug_alco_treat_tr_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.drug.alco.treat.tr<-lm(LRIRR ~ 1 + mods$reported_drug_or_alcohol_treatment_treat,weights = 1/v)

###Extract moderator Beta
mod_rep_drug_alco_treat_tr_LRIRR<-summary(model.lm.mod.rep.drug.alco.treat.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_drug_alco_treat_tr_se<-(summary(model.lm.mod.rep.drug.alco.treat.tr)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_drug_alco_treat_tr_res_df<- model.lm.mod.rep.drug.alco.treat.tr$df.residual

###Calculate Qres for the model
mod_rep_drug_alco_treat_tr_Y <- as.matrix(mods$LRIRR)
mod_rep_drug_alco_treat_tr_W <- diag(1/mods$RIRR_var)
mod_rep_drug_alco_treat_tr_X <- model.matrix(model.lm.mod.rep.drug.alco.treat.tr)
mod_rep_drug_alco_treat_tr_P <- mod_rep_drug_alco_treat_tr_W - mod_rep_drug_alco_treat_tr_W %*% mod_rep_drug_alco_treat_tr_X %*% solve(
  t(mod_rep_drug_alco_treat_tr_X) %*% mod_rep_drug_alco_treat_tr_W %*% mod_rep_drug_alco_treat_tr_X) %*% t(
    mod_rep_drug_alco_treat_tr_X) %*% mod_rep_drug_alco_treat_tr_W 

mod_rep_drug_alco_treat_tr_Qres <- max(0, c(crossprod(mod_rep_drug_alco_treat_tr_Y, mod_rep_drug_alco_treat_tr_P) %*% mod_rep_drug_alco_treat_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_drug_alco_treat_tr_list<-c(mod_rep_drug_alco_treat_tr_names,sg_rep_drug_alco_treat_tr_No_k, sg_rep_drug_alco_treat_tr_No_PLRIRR,sg_rep_drug_alco_treat_tr_No_PLRIRR_v, sg_rep_drug_alco_treat_tr_No_Q, sg_rep_drug_alco_treat_tr_No_Q_p, 
                                   sg_rep_drug_alco_treat_tr_No_I2,sg_rep_drug_alco_treat_tr_Yes_k, sg_rep_drug_alco_treat_tr_Yes_PLRIRR,sg_rep_drug_alco_treat_tr_Yes_PLRIRR_v, sg_rep_drug_alco_treat_tr_Yes_Q, 
                                   sg_rep_drug_alco_treat_tr_Yes_Q_p, sg_rep_drug_alco_treat_tr_Yes_I2, mod_rep_drug_alco_treat_tr_LRIRR, mod_rep_drug_alco_treat_tr_se,
                                   mod_rep_drug_alco_treat_tr_res_df, mod_rep_drug_alco_treat_tr_Qres)

names(mod_rep_drug_alco_treat_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                          "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                          "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_drug_alco_treat_tr_list)
########################################
### 8.5.6) Reported drug or alcohol treatment - Control
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_drug_alco_treat_con_names <-c('reported_drug_or_alcohol_treatment_cont','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_drug_alco_treat_con_No_es <-mods[mods$reported_drug_or_alcohol_treatment_cont== 'No', "LRIRR"]
sg_rep_drug_alco_treat_con_No_v <-mods[mods$reported_drug_or_alcohol_treatment_cont== 'No', "RIRR_var"]
sg_rep_drug_alco_treat_con_No_k <-length(sg_rep_drug_alco_treat_con_No_es)
sg_rep_drug_alco_treat_con_No_PLRIRR <-(sum(sg_rep_drug_alco_treat_con_No_es/sg_rep_drug_alco_treat_con_No_v))/(sum(1/sg_rep_drug_alco_treat_con_No_v))
sg_rep_drug_alco_treat_con_No_PLRIRR_v <- 1/sum(1/sg_rep_drug_alco_treat_con_No_v)
sg_rep_drug_alco_treat_con_No_Q <-sum((sg_rep_drug_alco_treat_con_No_es-sg_rep_drug_alco_treat_con_No_PLRIRR)^2/sg_rep_drug_alco_treat_con_No_v)
sg_rep_drug_alco_treat_con_No_Q_p <-pchisq(sg_rep_drug_alco_treat_con_No_Q, df=(length(sg_rep_drug_alco_treat_con_No_es)-1), lower.tail = FALSE)
sg_rep_drug_alco_treat_con_No_I2 <-((sg_rep_drug_alco_treat_con_No_Q-(length(sg_rep_drug_alco_treat_con_No_es)-1))/sg_rep_drug_alco_treat_con_No_Q)*100 

sg_rep_drug_alco_treat_con_Yes_es <-mods[mods$reported_drug_or_alcohol_treatment_cont== 'Yes', "LRIRR"]
sg_rep_drug_alco_treat_con_Yes_v <-mods[mods$reported_drug_or_alcohol_treatment_cont== 'Yes', "RIRR_var"]
sg_rep_drug_alco_treat_con_Yes_k <-length(sg_rep_drug_alco_treat_con_Yes_es)
sg_rep_drug_alco_treat_con_Yes_PLRIRR <-(sum(sg_rep_drug_alco_treat_con_Yes_es/sg_rep_drug_alco_treat_con_Yes_v))/(sum(1/sg_rep_drug_alco_treat_con_Yes_v))
sg_rep_drug_alco_treat_con_Yes_PLRIRR_v <- 1/sum(1/sg_rep_drug_alco_treat_con_Yes_v)
sg_rep_drug_alco_treat_con_Yes_Q <-sum((sg_rep_drug_alco_treat_con_Yes_es-sg_rep_drug_alco_treat_con_Yes_PLRIRR)^2/sg_rep_drug_alco_treat_con_Yes_v)
sg_rep_drug_alco_treat_con_Yes_Q_p <-pchisq(sg_rep_drug_alco_treat_con_Yes_Q, df=(length(sg_rep_drug_alco_treat_con_Yes_es)-1), lower.tail = FALSE)
sg_rep_drug_alco_treat_con_Yes_I2 <-((sg_rep_drug_alco_treat_con_Yes_Q-(length(sg_rep_drug_alco_treat_con_Yes_es)-1))/sg_rep_drug_alco_treat_con_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.drug.alco.treat.con<-lm(LRIRR ~ 1 + mods$reported_drug_or_alcohol_treatment_cont,weights = 1/v)

###Extract moderator Beta
mod_rep_drug_alco_treat_con_LRIRR<-summary(model.lm.mod.rep.drug.alco.treat.con)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_drug_alco_treat_con_se<-(summary(model.lm.mod.rep.drug.alco.treat.con)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_drug_alco_treat_con_res_df<- model.lm.mod.rep.drug.alco.treat.con$df.residual

###Calculate Qres for the model
mod_rep_drug_alco_treat_con_Y <- as.matrix(mods$LRIRR)
mod_rep_drug_alco_treat_con_W <- diag(1/mods$RIRR_var)
mod_rep_drug_alco_treat_con_X <- model.matrix(model.lm.mod.rep.drug.alco.treat.con)
mod_rep_drug_alco_treat_con_P <- mod_rep_drug_alco_treat_con_W - mod_rep_drug_alco_treat_con_W %*% mod_rep_drug_alco_treat_con_X %*% solve(
  t(mod_rep_drug_alco_treat_con_X) %*% mod_rep_drug_alco_treat_con_W %*% mod_rep_drug_alco_treat_con_X) %*% t(
    mod_rep_drug_alco_treat_con_X) %*% mod_rep_drug_alco_treat_con_W 

mod_rep_drug_alco_treat_con_Qres <- max(0, c(crossprod(mod_rep_drug_alco_treat_con_Y, mod_rep_drug_alco_treat_con_P) %*% mod_rep_drug_alco_treat_con_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_drug_alco_treat_con_list<-c(mod_rep_drug_alco_treat_con_names,sg_rep_drug_alco_treat_con_No_k, sg_rep_drug_alco_treat_con_No_PLRIRR,sg_rep_drug_alco_treat_con_No_PLRIRR_v, sg_rep_drug_alco_treat_con_No_Q, sg_rep_drug_alco_treat_con_No_Q_p, 
                                    sg_rep_drug_alco_treat_con_No_I2,sg_rep_drug_alco_treat_con_Yes_k, sg_rep_drug_alco_treat_con_Yes_PLRIRR,sg_rep_drug_alco_treat_con_Yes_PLRIRR_v, sg_rep_drug_alco_treat_con_Yes_Q, 
                                    sg_rep_drug_alco_treat_con_Yes_Q_p, sg_rep_drug_alco_treat_con_Yes_I2, mod_rep_drug_alco_treat_con_LRIRR, mod_rep_drug_alco_treat_con_se,
                                    mod_rep_drug_alco_treat_con_res_df, mod_rep_drug_alco_treat_con_Qres)

names(mod_rep_drug_alco_treat_con_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                           "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                           "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_drug_alco_treat_con_list)
########################################
### 8.5.7) Reported drug or alcohol testing - Treatment
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_drug_alco_test_tr_names <-c('reported_drug_or_alcohol_testing_treat','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_drug_alco_test_tr_No_es <-mods[mods$reported_drug_or_alcohol_testing_treat== 'No', "LRIRR"]
sg_rep_drug_alco_test_tr_No_v <-mods[mods$reported_drug_or_alcohol_testing_treat== 'No', "RIRR_var"]
sg_rep_drug_alco_test_tr_No_k <-length(sg_rep_drug_alco_test_tr_No_es)
sg_rep_drug_alco_test_tr_No_PLRIRR <-(sum(sg_rep_drug_alco_test_tr_No_es/sg_rep_drug_alco_test_tr_No_v))/(sum(1/sg_rep_drug_alco_test_tr_No_v))
sg_rep_drug_alco_test_tr_No_PLRIRR_v <- 1/sum(1/sg_rep_drug_alco_test_tr_No_v)
sg_rep_drug_alco_test_tr_No_Q <-sum((sg_rep_drug_alco_test_tr_No_es-sg_rep_drug_alco_test_tr_No_PLRIRR)^2/sg_rep_drug_alco_test_tr_No_v)
sg_rep_drug_alco_test_tr_No_Q_p <-pchisq(sg_rep_drug_alco_test_tr_No_Q, df=(length(sg_rep_drug_alco_test_tr_No_es)-1), lower.tail = FALSE)
sg_rep_drug_alco_test_tr_No_I2 <-((sg_rep_drug_alco_test_tr_No_Q-(length(sg_rep_drug_alco_test_tr_No_es)-1))/sg_rep_drug_alco_test_tr_No_Q)*100 

sg_rep_drug_alco_test_tr_Yes_es <-mods[mods$reported_drug_or_alcohol_testing_treat== 'Yes', "LRIRR"]
sg_rep_drug_alco_test_tr_Yes_v <-mods[mods$reported_drug_or_alcohol_testing_treat== 'Yes', "RIRR_var"]
sg_rep_drug_alco_test_tr_Yes_k <-length(sg_rep_drug_alco_test_tr_Yes_es)
sg_rep_drug_alco_test_tr_Yes_PLRIRR <-(sum(sg_rep_drug_alco_test_tr_Yes_es/sg_rep_drug_alco_test_tr_Yes_v))/(sum(1/sg_rep_drug_alco_test_tr_Yes_v))
sg_rep_drug_alco_test_tr_Yes_PLRIRR_v <- 1/sum(1/sg_rep_drug_alco_test_tr_Yes_v)
sg_rep_drug_alco_test_tr_Yes_Q <-sum((sg_rep_drug_alco_test_tr_Yes_es-sg_rep_drug_alco_test_tr_Yes_PLRIRR)^2/sg_rep_drug_alco_test_tr_Yes_v)
sg_rep_drug_alco_test_tr_Yes_Q_p <-pchisq(sg_rep_drug_alco_test_tr_Yes_Q, df=(length(sg_rep_drug_alco_test_tr_Yes_es)-1), lower.tail = FALSE)
sg_rep_drug_alco_test_tr_Yes_I2 <-((sg_rep_drug_alco_test_tr_Yes_Q-(length(sg_rep_drug_alco_test_tr_Yes_es)-1))/sg_rep_drug_alco_test_tr_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.drug.alco.test.tr<-lm(LRIRR ~ 1 + mods$reported_drug_or_alcohol_testing_treat,weights = 1/v)

###Extract moderator Beta
mod_rep_drug_alco_test_tr_LRIRR<-summary(model.lm.mod.rep.drug.alco.test.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_drug_alco_test_tr_se<-(summary(model.lm.mod.rep.drug.alco.test.tr)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_drug_alco_test_tr_res_df<- model.lm.mod.rep.drug.alco.test.tr$df.residual

###Calculate Qres for the model
mod_rep_drug_alco_test_tr_Y <- as.matrix(mods$LRIRR)
mod_rep_drug_alco_test_tr_W <- diag(1/mods$RIRR_var)
mod_rep_drug_alco_test_tr_X <- model.matrix(model.lm.mod.rep.drug.alco.test.tr)
mod_rep_drug_alco_test_tr_P <- mod_rep_drug_alco_test_tr_W - mod_rep_drug_alco_test_tr_W %*% mod_rep_drug_alco_test_tr_X %*% solve(
  t(mod_rep_drug_alco_test_tr_X) %*% mod_rep_drug_alco_test_tr_W %*% mod_rep_drug_alco_test_tr_X) %*% t(
    mod_rep_drug_alco_test_tr_X) %*% mod_rep_drug_alco_test_tr_W 

mod_rep_drug_alco_test_tr_Qres <- max(0, c(crossprod(mod_rep_drug_alco_test_tr_Y, mod_rep_drug_alco_test_tr_P) %*% mod_rep_drug_alco_test_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_drug_alco_test_tr_list<-c(mod_rep_drug_alco_test_tr_names,sg_rep_drug_alco_test_tr_No_k, sg_rep_drug_alco_test_tr_No_PLRIRR,sg_rep_drug_alco_test_tr_No_PLRIRR_v, sg_rep_drug_alco_test_tr_No_Q, sg_rep_drug_alco_test_tr_No_Q_p, 
                                  sg_rep_drug_alco_test_tr_No_I2,sg_rep_drug_alco_test_tr_Yes_k, sg_rep_drug_alco_test_tr_Yes_PLRIRR,sg_rep_drug_alco_test_tr_Yes_PLRIRR_v, sg_rep_drug_alco_test_tr_Yes_Q, 
                                  sg_rep_drug_alco_test_tr_Yes_Q_p, sg_rep_drug_alco_test_tr_Yes_I2, mod_rep_drug_alco_test_tr_LRIRR, mod_rep_drug_alco_test_tr_se,
                                  mod_rep_drug_alco_test_tr_res_df, mod_rep_drug_alco_test_tr_Qres)

names(mod_rep_drug_alco_test_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                         "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                         "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_drug_alco_test_tr_list)
########################################
### 8.5.8) Reported drug or alcohol testing - Control
########################################
### Note: there was insufficient variance in this moderator to perform an analysis 
########################################
### 8.5.9) Reported psychological treatment - Treatment
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_psych_treat_tr_names <-c('reported_psych_treatment_treat','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_psych_treat_tr_No_es <-mods[mods$reported_psychiatric_treatment_treat== 'No', "LRIRR"]
sg_rep_psych_treat_tr_No_v <-mods[mods$reported_psychiatric_treatment_treat== 'No', "RIRR_var"]
sg_rep_psych_treat_tr_No_k <-length(sg_rep_psych_treat_tr_No_es)
sg_rep_psych_treat_tr_No_PLRIRR <-(sum(sg_rep_psych_treat_tr_No_es/sg_rep_psych_treat_tr_No_v))/(sum(1/sg_rep_psych_treat_tr_No_v))
sg_rep_psych_treat_tr_No_PLRIRR_v <- 1/sum(1/sg_rep_psych_treat_tr_No_v)
sg_rep_psych_treat_tr_No_Q <-sum((sg_rep_psych_treat_tr_No_es-sg_rep_psych_treat_tr_No_PLRIRR)^2/sg_rep_psych_treat_tr_No_v)
sg_rep_psych_treat_tr_No_Q_p <-pchisq(sg_rep_psych_treat_tr_No_Q, df=(length(sg_rep_psych_treat_tr_No_es)-1), lower.tail = FALSE)
sg_rep_psych_treat_tr_No_I2 <-((sg_rep_psych_treat_tr_No_Q-(length(sg_rep_psych_treat_tr_No_es)-1))/sg_rep_psych_treat_tr_No_Q)*100 

sg_rep_psych_treat_tr_Yes_es <-mods[mods$reported_psychiatric_treatment_treat== 'Yes', "LRIRR"]
sg_rep_psych_treat_tr_Yes_v <-mods[mods$reported_psychiatric_treatment_treat== 'Yes', "RIRR_var"]
sg_rep_psych_treat_tr_Yes_k <-length(sg_rep_psych_treat_tr_Yes_es)
sg_rep_psych_treat_tr_Yes_PLRIRR <-(sum(sg_rep_psych_treat_tr_Yes_es/sg_rep_psych_treat_tr_Yes_v))/(sum(1/sg_rep_psych_treat_tr_Yes_v))
sg_rep_psych_treat_tr_Yes_PLRIRR_v <- 1/sum(1/sg_rep_psych_treat_tr_Yes_v)
sg_rep_psych_treat_tr_Yes_Q <-sum((sg_rep_psych_treat_tr_Yes_es-sg_rep_psych_treat_tr_Yes_PLRIRR)^2/sg_rep_psych_treat_tr_Yes_v)
sg_rep_psych_treat_tr_Yes_Q_p <-pchisq(sg_rep_psych_treat_tr_Yes_Q, df=(length(sg_rep_psych_treat_tr_Yes_es)-1), lower.tail = FALSE)
sg_rep_psych_treat_tr_Yes_I2 <-((sg_rep_psych_treat_tr_Yes_Q-(length(sg_rep_psych_treat_tr_Yes_es)-1))/sg_rep_psych_treat_tr_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.psych.treat.tr<-lm(LRIRR ~ 1 + mods$reported_psychiatric_treatment_treat,weights = 1/v)

###Extract moderator Beta
mod_rep_psych_treat_tr_LRIRR<-summary(model.lm.mod.rep.psych.treat.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_psych_treat_tr_se<-(summary(model.lm.mod.rep.psych.treat.tr)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_psych_treat_tr_res_df<- model.lm.mod.rep.psych.treat.tr$df.residual

###Calculate Qres for the model
mod_rep_psych_treat_tr_Y <- as.matrix(mods$LRIRR)
mod_rep_psych_treat_tr_W <- diag(1/mods$RIRR_var)
mod_rep_psych_treat_tr_X <- model.matrix(model.lm.mod.rep.psych.treat.tr)
mod_rep_psych_treat_tr_P <- mod_rep_psych_treat_tr_W - mod_rep_psych_treat_tr_W %*% mod_rep_psych_treat_tr_X %*% solve(
  t(mod_rep_psych_treat_tr_X) %*% mod_rep_psych_treat_tr_W %*% mod_rep_psych_treat_tr_X) %*% t(
    mod_rep_psych_treat_tr_X) %*% mod_rep_psych_treat_tr_W 

mod_rep_psych_treat_tr_Qres <- max(0, c(crossprod(mod_rep_psych_treat_tr_Y, mod_rep_psych_treat_tr_P) %*% mod_rep_psych_treat_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_psych_treat_tr_list<-c(mod_rep_psych_treat_tr_names,sg_rep_psych_treat_tr_No_k, sg_rep_psych_treat_tr_No_PLRIRR,sg_rep_psych_treat_tr_No_PLRIRR_v, sg_rep_psych_treat_tr_No_Q, sg_rep_psych_treat_tr_No_Q_p, 
                               sg_rep_psych_treat_tr_No_I2,sg_rep_psych_treat_tr_Yes_k, sg_rep_psych_treat_tr_Yes_PLRIRR,sg_rep_psych_treat_tr_Yes_PLRIRR_v, sg_rep_psych_treat_tr_Yes_Q, 
                               sg_rep_psych_treat_tr_Yes_Q_p, sg_rep_psych_treat_tr_Yes_I2, mod_rep_psych_treat_tr_LRIRR, mod_rep_psych_treat_tr_se,
                               mod_rep_psych_treat_tr_res_df, mod_rep_psych_treat_tr_Qres)

names(mod_rep_psych_treat_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                      "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                      "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_psych_treat_tr_list)

########################################
### 8.5.10) Reported psychological treatment - Control
########################################
### Note: there was insufficient variance in this moderator to perform an analysis 
########################################
### 8.5.11) Reported counselling - Treatment
########################################
# Note: this variable separates the subgroups in the same way as '8.5.2) Reported measuring risk - Control' and so will be excluded from the mean moderator analysis.
all.equal(table(mods$reported_counseling_treat, mods$studyid), table(mods$reported_measuring_risk_cont, mods$studyid))
########################################
### 8.5.12) Reported counselling - Control
########################################
### Note: there was insufficient variance in this moderator to perform an analysis 
########################################
### 8.5.13) Reported support groups - Treatment
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_supp_groups_tr_names <-c('reported_support_groups_treat','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_supp_groups_tr_No_es <-mods[mods$reported_support_groups_treat== 'No', "LRIRR"]
sg_rep_supp_groups_tr_No_v <-mods[mods$reported_support_groups_treat== 'No', "RIRR_var"]
sg_rep_supp_groups_tr_No_k <-length(sg_rep_supp_groups_tr_No_es)
sg_rep_supp_groups_tr_No_PLRIRR <-(sum(sg_rep_supp_groups_tr_No_es/sg_rep_supp_groups_tr_No_v))/(sum(1/sg_rep_supp_groups_tr_No_v))
sg_rep_supp_groups_tr_No_PLRIRR_v <- 1/sum(1/sg_rep_supp_groups_tr_No_v)
sg_rep_supp_groups_tr_No_Q <-sum((sg_rep_supp_groups_tr_No_es-sg_rep_supp_groups_tr_No_PLRIRR)^2/sg_rep_supp_groups_tr_No_v)
sg_rep_supp_groups_tr_No_Q_p <-pchisq(sg_rep_supp_groups_tr_No_Q, df=(length(sg_rep_supp_groups_tr_No_es)-1), lower.tail = FALSE)
sg_rep_supp_groups_tr_No_I2 <-((sg_rep_supp_groups_tr_No_Q-(length(sg_rep_supp_groups_tr_No_es)-1))/sg_rep_supp_groups_tr_No_Q)*100 

sg_rep_supp_groups_tr_Yes_es <-mods[mods$reported_support_groups_treat== 'Yes', "LRIRR"]
sg_rep_supp_groups_tr_Yes_v <-mods[mods$reported_support_groups_treat== 'Yes', "RIRR_var"]
sg_rep_supp_groups_tr_Yes_k <-length(sg_rep_supp_groups_tr_Yes_es)
sg_rep_supp_groups_tr_Yes_PLRIRR <-(sum(sg_rep_supp_groups_tr_Yes_es/sg_rep_supp_groups_tr_Yes_v))/(sum(1/sg_rep_supp_groups_tr_Yes_v))
sg_rep_supp_groups_tr_Yes_PLRIRR_v <- 1/sum(1/sg_rep_supp_groups_tr_Yes_v)
sg_rep_supp_groups_tr_Yes_Q <-sum((sg_rep_supp_groups_tr_Yes_es-sg_rep_supp_groups_tr_Yes_PLRIRR)^2/sg_rep_supp_groups_tr_Yes_v)
sg_rep_supp_groups_tr_Yes_Q_p <-pchisq(sg_rep_supp_groups_tr_Yes_Q, df=(length(sg_rep_supp_groups_tr_Yes_es)-1), lower.tail = FALSE)
sg_rep_supp_groups_tr_Yes_I2 <-((sg_rep_supp_groups_tr_Yes_Q-(length(sg_rep_supp_groups_tr_Yes_es)-1))/sg_rep_supp_groups_tr_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.supp.groups.tr<-lm(LRIRR ~ 1 + mods$reported_support_groups_treat,weights = 1/v)

###Extract moderator Beta
mod_rep_supp_groups_tr_LRIRR<-summary(model.lm.mod.rep.supp.groups.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_supp_groups_tr_se<-(summary(model.lm.mod.rep.supp.groups.tr)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_supp_groups_tr_res_df<- model.lm.mod.rep.supp.groups.tr$df.residual

###Calculate Qres for the model
mod_rep_supp_groups_tr_Y <- as.matrix(mods$LRIRR)
mod_rep_supp_groups_tr_W <- diag(1/mods$RIRR_var)
mod_rep_supp_groups_tr_X <- model.matrix(model.lm.mod.rep.supp.groups.tr)
mod_rep_supp_groups_tr_P <- mod_rep_supp_groups_tr_W - mod_rep_supp_groups_tr_W %*% mod_rep_supp_groups_tr_X %*% solve(
  t(mod_rep_supp_groups_tr_X) %*% mod_rep_supp_groups_tr_W %*% mod_rep_supp_groups_tr_X) %*% t(
    mod_rep_supp_groups_tr_X) %*% mod_rep_supp_groups_tr_W 

mod_rep_supp_groups_tr_Qres <- max(0, c(crossprod(mod_rep_supp_groups_tr_Y, mod_rep_supp_groups_tr_P) %*% mod_rep_supp_groups_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_supp_groups_tr_list<-c(mod_rep_supp_groups_tr_names,sg_rep_supp_groups_tr_No_k, sg_rep_supp_groups_tr_No_PLRIRR,sg_rep_supp_groups_tr_No_PLRIRR_v, sg_rep_supp_groups_tr_No_Q, sg_rep_supp_groups_tr_No_Q_p, 
                               sg_rep_supp_groups_tr_No_I2,sg_rep_supp_groups_tr_Yes_k, sg_rep_supp_groups_tr_Yes_PLRIRR,sg_rep_supp_groups_tr_Yes_PLRIRR_v, sg_rep_supp_groups_tr_Yes_Q, 
                               sg_rep_supp_groups_tr_Yes_Q_p, sg_rep_supp_groups_tr_Yes_I2, mod_rep_supp_groups_tr_LRIRR, mod_rep_supp_groups_tr_se,
                               mod_rep_supp_groups_tr_res_df, mod_rep_supp_groups_tr_Qres)

names(mod_rep_supp_groups_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                      "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                      "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_supp_groups_tr_list)
########################################
### 8.5.14) Reported support groups - Control
########################################
### Note: there was insufficient variance in this moderator to perform an analysis 
########################################
### 8.5.15) Reported probation - Control
########################################
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_probation_con_names <-c('reported_probation_cont','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_probation_con_No_es <-mods[mods$reported_probation_cont== 'No', "LRIRR"]
sg_rep_probation_con_No_v <-mods[mods$reported_probation_cont== 'No', "RIRR_var"]
sg_rep_probation_con_No_k <-length(sg_rep_probation_con_No_es)
sg_rep_probation_con_No_PLRIRR <-(sum(sg_rep_probation_con_No_es/sg_rep_probation_con_No_v))/(sum(1/sg_rep_probation_con_No_v))
sg_rep_probation_con_No_PLRIRR_v <- 1/sum(1/sg_rep_probation_con_No_v)
sg_rep_probation_con_No_Q <-sum((sg_rep_probation_con_No_es-sg_rep_probation_con_No_PLRIRR)^2/sg_rep_probation_con_No_v)
sg_rep_probation_con_No_Q_p <-pchisq(sg_rep_probation_con_No_Q, df=(length(sg_rep_probation_con_No_es)-1), lower.tail = FALSE)
sg_rep_probation_con_No_I2 <-((sg_rep_probation_con_No_Q-(length(sg_rep_probation_con_No_es)-1))/sg_rep_probation_con_No_Q)*100 

sg_rep_probation_con_Yes_es <-mods[mods$reported_probation_cont== 'Yes', "LRIRR"]
sg_rep_probation_con_Yes_v <-mods[mods$reported_probation_cont== 'Yes', "RIRR_var"]
sg_rep_probation_con_Yes_k <-length(sg_rep_probation_con_Yes_es)
sg_rep_probation_con_Yes_PLRIRR <-(sum(sg_rep_probation_con_Yes_es/sg_rep_probation_con_Yes_v))/(sum(1/sg_rep_probation_con_Yes_v))
sg_rep_probation_con_Yes_PLRIRR_v <- 1/sum(1/sg_rep_probation_con_Yes_v)
sg_rep_probation_con_Yes_Q <-sum((sg_rep_probation_con_Yes_es-sg_rep_probation_con_Yes_PLRIRR)^2/sg_rep_probation_con_Yes_v)
sg_rep_probation_con_Yes_Q_p <-pchisq(sg_rep_probation_con_Yes_Q, df=(length(sg_rep_probation_con_Yes_es)-1), lower.tail = FALSE)
sg_rep_probation_con_Yes_I2 <-((sg_rep_probation_con_Yes_Q-(length(sg_rep_probation_con_Yes_es)-1))/sg_rep_probation_con_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.probation.con<-lm(LRIRR ~ 1 + mods$reported_probation_cont,weights = 1/v)

###Extract moderator Beta
mod_rep_probation_con_LRIRR<-summary(model.lm.mod.rep.probation.con)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_probation_con_se<-(summary(model.lm.mod.rep.probation.con)$coefficients[2, 2])

###Create an object for the residual degrees of freedom. 
mod_rep_probation_con_res_df<- model.lm.mod.rep.probation.con$df.residual

###Calculate Qres for the model
mod_rep_probation_con_Y <- as.matrix(mods$LRIRR)
mod_rep_probation_con_W <- diag(1/mods$RIRR_var)
mod_rep_probation_con_X <- model.matrix(model.lm.mod.rep.probation.con)
mod_rep_probation_con_P <- mod_rep_probation_con_W - mod_rep_probation_con_W %*% mod_rep_probation_con_X %*% solve(
  t(mod_rep_probation_con_X) %*% mod_rep_probation_con_W %*% mod_rep_probation_con_X) %*% t(
    mod_rep_probation_con_X) %*% mod_rep_probation_con_W 

mod_rep_probation_con_Qres <- max(0, c(crossprod(mod_rep_probation_con_Y, mod_rep_probation_con_P) %*% mod_rep_probation_con_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_probation_con_list<-c(mod_rep_probation_con_names,sg_rep_probation_con_No_k, sg_rep_probation_con_No_PLRIRR,sg_rep_probation_con_No_PLRIRR_v, sg_rep_probation_con_No_Q, sg_rep_probation_con_No_Q_p, 
                              sg_rep_probation_con_No_I2,sg_rep_probation_con_Yes_k, sg_rep_probation_con_Yes_PLRIRR,sg_rep_probation_con_Yes_PLRIRR_v, sg_rep_probation_con_Yes_Q, 
                              sg_rep_probation_con_Yes_Q_p, sg_rep_probation_con_Yes_I2, mod_rep_probation_con_LRIRR, mod_rep_probation_con_se,
                              mod_rep_probation_con_res_df, mod_rep_probation_con_Qres)

names(mod_rep_probation_con_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                     "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                     "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_probation_con_list)
########################################
### 8.5.16) Reported other treatment - Aftercare - Treatment
########################################
###Create a variable that seperates studies that reported aftercare from studies that did not. 
mods$other_treatment_treat_aftercare<- mods$other_treatment_treat
mods$other_treatment_treat_aftercare<- mods$other_treatment_treat_aftercare == 'aftercare'
mods$other_treatment_treat_aftercare<- na_if(mods$other_treatment_treat_aftercare, 'FALSE')
mods$other_treatment_treat_aftercare[is.na(mods$other_treatment_treat_aftercare)] <- 'FALSE'
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_aftercare_tr_names <-c('reported_aftercare_treat','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_aftercare_tr_No_es <-mods[mods$other_treatment_treat_aftercare== 'FALSE', "LRIRR"]
sg_rep_aftercare_tr_No_v <-mods[mods$other_treatment_treat_aftercare== 'FALSE', "RIRR_var"]
sg_rep_aftercare_tr_No_k <-length(sg_rep_aftercare_tr_No_es)
sg_rep_aftercare_tr_No_PLRIRR <-(sum(sg_rep_aftercare_tr_No_es/sg_rep_aftercare_tr_No_v))/(sum(1/sg_rep_aftercare_tr_No_v))
sg_rep_aftercare_tr_No_PLRIRR_v <- 1/sum(1/sg_rep_aftercare_tr_No_v)
sg_rep_aftercare_tr_No_Q <-sum((sg_rep_aftercare_tr_No_es-sg_rep_aftercare_tr_No_PLRIRR)^2/sg_rep_aftercare_tr_No_v)
sg_rep_aftercare_tr_No_Q_p <-pchisq(sg_rep_aftercare_tr_No_Q, df=(length(sg_rep_aftercare_tr_No_es)-1), lower.tail = FALSE)
sg_rep_aftercare_tr_No_I2 <-((sg_rep_aftercare_tr_No_Q-(length(sg_rep_aftercare_tr_No_es)-1))/sg_rep_aftercare_tr_No_Q)*100 

sg_rep_aftercare_tr_Yes_es <-mods[mods$other_treatment_treat_aftercare== 'TRUE', "LRIRR"]
sg_rep_aftercare_tr_Yes_v <-mods[mods$other_treatment_treat_aftercare== 'TRUE', "RIRR_var"]
sg_rep_aftercare_tr_Yes_k <-length(sg_rep_aftercare_tr_Yes_es)
sg_rep_aftercare_tr_Yes_PLRIRR <-(sum(sg_rep_aftercare_tr_Yes_es/sg_rep_aftercare_tr_Yes_v))/(sum(1/sg_rep_aftercare_tr_Yes_v))
sg_rep_aftercare_tr_Yes_PLRIRR_v <- 1/sum(1/sg_rep_aftercare_tr_Yes_v)
sg_rep_aftercare_tr_Yes_Q <-sum((sg_rep_aftercare_tr_Yes_es-sg_rep_aftercare_tr_Yes_PLRIRR)^2/sg_rep_aftercare_tr_Yes_v)
sg_rep_aftercare_tr_Yes_Q_p <-pchisq(sg_rep_aftercare_tr_Yes_Q, df=(length(sg_rep_aftercare_tr_Yes_es)-1), lower.tail = FALSE)
sg_rep_aftercare_tr_Yes_I2 <-((sg_rep_aftercare_tr_Yes_Q-(length(sg_rep_aftercare_tr_Yes_es)-1))/sg_rep_aftercare_tr_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.aftercare.tr<-lm(LRIRR ~ 1 + mods$other_treatment_treat_aftercare,weights = 1/v)

###Extract moderator Beta
mod_rep_aftercare_tr_LRIRR<-summary(model.lm.mod.rep.aftercare.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_aftercare_tr_se<-summary(model.lm.mod.rep.aftercare.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_rep_aftercare_tr_res_df<- model.lm.mod.rep.aftercare.tr$df.residual

###Calculate Qres for the model
mod_rep_aftercare_tr_Y <- as.matrix(mods$LRIRR)
mod_rep_aftercare_tr_W <- diag(1/mods$RIRR_var)
mod_rep_aftercare_tr_X <- model.matrix(model.lm.mod.rep.aftercare.tr)
mod_rep_aftercare_tr_P <- mod_rep_aftercare_tr_W - mod_rep_aftercare_tr_W %*% mod_rep_aftercare_tr_X %*% solve(
  t(mod_rep_aftercare_tr_X) %*% mod_rep_aftercare_tr_W %*% mod_rep_aftercare_tr_X) %*% t(
    mod_rep_aftercare_tr_X) %*% mod_rep_aftercare_tr_W 

mod_rep_aftercare_tr_Qres <- max(0, c(crossprod(mod_rep_aftercare_tr_Y, mod_rep_aftercare_tr_P) %*% mod_rep_aftercare_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_aftercare_tr_list<-c(mod_rep_aftercare_tr_names,sg_rep_aftercare_tr_No_k, sg_rep_aftercare_tr_No_PLRIRR,sg_rep_aftercare_tr_No_PLRIRR_v, sg_rep_aftercare_tr_No_Q, sg_rep_aftercare_tr_No_Q_p, 
                             sg_rep_aftercare_tr_No_I2,sg_rep_aftercare_tr_Yes_k, sg_rep_aftercare_tr_Yes_PLRIRR,sg_rep_aftercare_tr_Yes_PLRIRR_v, sg_rep_aftercare_tr_Yes_Q, 
                             sg_rep_aftercare_tr_Yes_Q_p, sg_rep_aftercare_tr_Yes_I2, mod_rep_aftercare_tr_LRIRR, mod_rep_aftercare_tr_se,
                             mod_rep_aftercare_tr_res_df, mod_rep_aftercare_tr_Qres)

names(mod_rep_aftercare_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                    "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                    "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_aftercare_tr_list)
########################################
### 8.5.17) Reported other treatment - Relapse prevention - Treatment
########################################
###Create a variable that seperates studies that reported relapse prevention from studies that did not. 
mods$other_treatment_treat_relap_prev<- mods$other_treatment_treat
mods$other_treatment_treat_relap_prev<- mods$other_treatment_treat_relap_prev == 'relapse prevention'
mods$other_treatment_treat_relap_prev<- na_if(mods$other_treatment_treat_relap_prev, 'FALSE')
mods$other_treatment_treat_relap_prev[is.na(mods$other_treatment_treat_relap_prev)] <- 'FALSE'
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_relap_prev_tr_names <-c('reported_relapse_prevention_treat','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_relap_prev_tr_No_es <-mods[mods$other_treatment_treat_relap_prev== 'FALSE', "LRIRR"]
sg_rep_relap_prev_tr_No_v <-mods[mods$other_treatment_treat_relap_prev== 'FALSE', "RIRR_var"]
sg_rep_relap_prev_tr_No_k <-length(sg_rep_relap_prev_tr_No_es)
sg_rep_relap_prev_tr_No_PLRIRR <-(sum(sg_rep_relap_prev_tr_No_es/sg_rep_relap_prev_tr_No_v))/(sum(1/sg_rep_relap_prev_tr_No_v))
sg_rep_relap_prev_tr_No_PLRIRR_v <- 1/sum(1/sg_rep_relap_prev_tr_No_v)
sg_rep_relap_prev_tr_No_Q <-sum((sg_rep_relap_prev_tr_No_es-sg_rep_relap_prev_tr_No_PLRIRR)^2/sg_rep_relap_prev_tr_No_v)
sg_rep_relap_prev_tr_No_Q_p <-pchisq(sg_rep_relap_prev_tr_No_Q, df=(length(sg_rep_relap_prev_tr_No_es)-1), lower.tail = FALSE)
sg_rep_relap_prev_tr_No_I2 <-((sg_rep_relap_prev_tr_No_Q-(length(sg_rep_relap_prev_tr_No_es)-1))/sg_rep_relap_prev_tr_No_Q)*100 

sg_rep_relap_prev_tr_Yes_es <-mods[mods$other_treatment_treat_relap_prev== 'TRUE', "LRIRR"]
sg_rep_relap_prev_tr_Yes_v <-mods[mods$other_treatment_treat_relap_prev== 'TRUE', "RIRR_var"]
sg_rep_relap_prev_tr_Yes_k <-length(sg_rep_relap_prev_tr_Yes_es)
sg_rep_relap_prev_tr_Yes_PLRIRR <-(sum(sg_rep_relap_prev_tr_Yes_es/sg_rep_relap_prev_tr_Yes_v))/(sum(1/sg_rep_relap_prev_tr_Yes_v))
sg_rep_relap_prev_tr_Yes_PLRIRR_v <- 1/sum(1/sg_rep_relap_prev_tr_Yes_v)
sg_rep_relap_prev_tr_Yes_Q <-sum((sg_rep_relap_prev_tr_Yes_es-sg_rep_relap_prev_tr_Yes_PLRIRR)^2/sg_rep_relap_prev_tr_Yes_v)
sg_rep_relap_prev_tr_Yes_Q_p <-pchisq(sg_rep_relap_prev_tr_Yes_Q, df=(length(sg_rep_relap_prev_tr_Yes_es)-1), lower.tail = FALSE)
sg_rep_relap_prev_tr_Yes_I2 <-((sg_rep_relap_prev_tr_Yes_Q-(length(sg_rep_relap_prev_tr_Yes_es)-1))/sg_rep_relap_prev_tr_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.relap.prev.tr<-lm(LRIRR ~ 1 + mods$other_treatment_treat_relap_prev,weights = 1/v)

###Extract moderator Beta
mod_rep_relap_prev_tr_LRIRR<-summary(model.lm.mod.rep.relap.prev.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_relap_prev_tr_se<-summary(model.lm.mod.rep.relap.prev.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_rep_relap_prev_tr_res_df<- model.lm.mod.rep.relap.prev.tr$df.residual

###Calculate Qres for the model
mod_rep_relap_prev_tr_Y <- as.matrix(mods$LRIRR)
mod_rep_relap_prev_tr_W <- diag(1/mods$RIRR_var)
mod_rep_relap_prev_tr_X <- model.matrix(model.lm.mod.rep.relap.prev.tr)
mod_rep_relap_prev_tr_P <- mod_rep_relap_prev_tr_W - mod_rep_relap_prev_tr_W %*% mod_rep_relap_prev_tr_X %*% solve(
  t(mod_rep_relap_prev_tr_X) %*% mod_rep_relap_prev_tr_W %*% mod_rep_relap_prev_tr_X) %*% t(
    mod_rep_relap_prev_tr_X) %*% mod_rep_relap_prev_tr_W 

mod_rep_relap_prev_tr_Qres <- max(0, c(crossprod(mod_rep_relap_prev_tr_Y, mod_rep_relap_prev_tr_P) %*% mod_rep_relap_prev_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_relap_prev_tr_list<-c(mod_rep_relap_prev_tr_names,sg_rep_relap_prev_tr_No_k, sg_rep_relap_prev_tr_No_PLRIRR,sg_rep_relap_prev_tr_No_PLRIRR_v, sg_rep_relap_prev_tr_No_Q, sg_rep_relap_prev_tr_No_Q_p, 
                              sg_rep_relap_prev_tr_No_I2,sg_rep_relap_prev_tr_Yes_k, sg_rep_relap_prev_tr_Yes_PLRIRR,sg_rep_relap_prev_tr_Yes_PLRIRR_v, sg_rep_relap_prev_tr_Yes_Q, 
                              sg_rep_relap_prev_tr_Yes_Q_p, sg_rep_relap_prev_tr_Yes_I2, mod_rep_relap_prev_tr_LRIRR, mod_rep_relap_prev_tr_se,
                              mod_rep_relap_prev_tr_res_df, mod_rep_relap_prev_tr_Qres)

names(mod_rep_relap_prev_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                     "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                     "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_relap_prev_tr_list)
########################################
### 8.5.18) Reported other treatment - Self-help - Treatment
########################################
###Create a variable that seperates studies that reported 'relapse prevention'self-help' from studies that did not. 
mods$other_treatment_treat_self_help<- mods$other_treatment_treat
mods$other_treatment_treat_self_help<- mods$other_treatment_treat_self_help == 'self-help'
mods$other_treatment_treat_self_help<- na_if(mods$other_treatment_treat_self_help, 'FALSE')
mods$other_treatment_treat_self_help[is.na(mods$other_treatment_treat_self_help)] <- 'FALSE'
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_self_help_tr_names <-c('reported_self_help_treat','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_self_help_tr_No_es <-mods[mods$other_treatment_treat_self_help== 'FALSE', "LRIRR"]
sg_rep_self_help_tr_No_v <-mods[mods$other_treatment_treat_self_help== 'FALSE', "RIRR_var"]
sg_rep_self_help_tr_No_k <-length(sg_rep_self_help_tr_No_es)
sg_rep_self_help_tr_No_PLRIRR <-(sum(sg_rep_self_help_tr_No_es/sg_rep_self_help_tr_No_v))/(sum(1/sg_rep_self_help_tr_No_v))
sg_rep_self_help_tr_No_PLRIRR_v <- 1/sum(1/sg_rep_self_help_tr_No_v)
sg_rep_self_help_tr_No_Q <-sum((sg_rep_self_help_tr_No_es-sg_rep_self_help_tr_No_PLRIRR)^2/sg_rep_self_help_tr_No_v)
sg_rep_self_help_tr_No_Q_p <-pchisq(sg_rep_self_help_tr_No_Q, df=(length(sg_rep_self_help_tr_No_es)-1), lower.tail = FALSE)
sg_rep_self_help_tr_No_I2 <-((sg_rep_self_help_tr_No_Q-(length(sg_rep_self_help_tr_No_es)-1))/sg_rep_self_help_tr_No_Q)*100 

sg_rep_self_help_tr_Yes_es <-mods[mods$other_treatment_treat_self_help== 'TRUE', "LRIRR"]
sg_rep_self_help_tr_Yes_v <-mods[mods$other_treatment_treat_self_help== 'TRUE', "RIRR_var"]
sg_rep_self_help_tr_Yes_k <-length(sg_rep_self_help_tr_Yes_es)
sg_rep_self_help_tr_Yes_PLRIRR <-(sum(sg_rep_self_help_tr_Yes_es/sg_rep_self_help_tr_Yes_v))/(sum(1/sg_rep_self_help_tr_Yes_v))
sg_rep_self_help_tr_Yes_PLRIRR_v <- 1/sum(1/sg_rep_self_help_tr_Yes_v)
sg_rep_self_help_tr_Yes_Q <-sum((sg_rep_self_help_tr_Yes_es-sg_rep_self_help_tr_Yes_PLRIRR)^2/sg_rep_self_help_tr_Yes_v)
sg_rep_self_help_tr_Yes_Q_p <-pchisq(sg_rep_self_help_tr_Yes_Q, df=(length(sg_rep_self_help_tr_Yes_es)-1), lower.tail = FALSE)
sg_rep_self_help_tr_Yes_I2 <-((sg_rep_self_help_tr_Yes_Q-(length(sg_rep_self_help_tr_Yes_es)-1))/sg_rep_self_help_tr_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.self.help.tr<-lm(LRIRR ~ 1 + mods$other_treatment_treat_self_help,weights = 1/v)

###Extract moderator Beta
mod_rep_self_help_tr_LRIRR<-summary(model.lm.mod.rep.self.help.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_self_help_tr_se<-summary(model.lm.mod.rep.self.help.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_rep_self_help_tr_res_df<- model.lm.mod.rep.self.help.tr$df.residual

###Calculate Qres for the model
mod_rep_self_help_tr_Y <- as.matrix(mods$LRIRR)
mod_rep_self_help_tr_W <- diag(1/mods$RIRR_var)
mod_rep_self_help_tr_X <- model.matrix(model.lm.mod.rep.self.help.tr)
mod_rep_self_help_tr_P <- mod_rep_self_help_tr_W - mod_rep_self_help_tr_W %*% mod_rep_self_help_tr_X %*% solve(
  t(mod_rep_self_help_tr_X) %*% mod_rep_self_help_tr_W %*% mod_rep_self_help_tr_X) %*% t(
    mod_rep_self_help_tr_X) %*% mod_rep_self_help_tr_W 

mod_rep_self_help_tr_Qres <- max(0, c(crossprod(mod_rep_self_help_tr_Y, mod_rep_self_help_tr_P) %*% mod_rep_self_help_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_self_help_tr_list<-c(mod_rep_self_help_tr_names,sg_rep_self_help_tr_No_k, sg_rep_self_help_tr_No_PLRIRR,sg_rep_self_help_tr_No_PLRIRR_v, sg_rep_self_help_tr_No_Q, sg_rep_self_help_tr_No_Q_p, 
                             sg_rep_self_help_tr_No_I2,sg_rep_self_help_tr_Yes_k, sg_rep_self_help_tr_Yes_PLRIRR,sg_rep_self_help_tr_Yes_PLRIRR_v, sg_rep_self_help_tr_Yes_Q, 
                             sg_rep_self_help_tr_Yes_Q_p, sg_rep_self_help_tr_Yes_I2, mod_rep_self_help_tr_LRIRR, mod_rep_self_help_tr_se,
                             mod_rep_self_help_tr_res_df, mod_rep_self_help_tr_Qres)

names(mod_rep_self_help_tr_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                    "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                    "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_self_help_tr_list)
########################################
### 8.5.19) Total treatments required - Treatment
########################################
### Create an object for the row name, level_ref, and level_alternative names
tot_treat_req_tr_names<- c('total_treatments_required_treat', 'treatment', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                           'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                           'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.tot.treat.req.tr<-lm(LRIRR ~ 1 + mods$total_treatments_required_treat,weights = 1/v)

###Extract moderator Beta
mod_tot_treat_req_tr_B_LRIRR<-summary(model.lm.mod.tot.treat.req.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_tot_treat_req_tr_B_se<-summary(model.lm.mod.tot.treat.req.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_tot_treat_req_tr_res_df<- model.lm.mod.tot.treat.req.tr$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$total_treatments_required_tr_is_na<- is.na(mods$total_treatments_required_treat)
mod_tot_treat_req_tr_es <- mods[mods$total_treatments_required_tr_is_na == 'FALSE', "LRIRR"]
mod_tot_treat_req_tr_v <- mods[mods$total_treatments_required_tr_is_na == 'FALSE', "RIRR_var"]

mod_tot_treat_req_tr_Y <- as.matrix(mod_tot_treat_req_tr_es)
mod_tot_treat_req_tr_W <- diag(1/mod_tot_treat_req_tr_v)
mod_tot_treat_req_tr_X <- model.matrix(model.lm.mod.tot.treat.req.tr)
mod_tot_treat_req_tr_P <- mod_tot_treat_req_tr_W - mod_tot_treat_req_tr_W %*% mod_tot_treat_req_tr_X %*% solve(
  t(mod_tot_treat_req_tr_X) %*% mod_tot_treat_req_tr_W %*% mod_tot_treat_req_tr_X) %*% t(
    mod_tot_treat_req_tr_X) %*% mod_tot_treat_req_tr_W 

mod_tot_treat_req_tr_Qres <- max(0, c(crossprod(mod_tot_treat_req_tr_Y, mod_tot_treat_req_tr_P) %*% mod_tot_treat_req_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
tot_treat_req_tr_names_list<- c(tot_treat_req_tr_names, mod_tot_treat_req_tr_B_LRIRR, mod_tot_treat_req_tr_B_se, 
                                mod_tot_treat_req_tr_res_df, mod_tot_treat_req_tr_Qres)

names(tot_treat_req_tr_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                      "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                      "mod_Qres")
moderator_output<- bind_rows(moderator_output, tot_treat_req_tr_names_list)
########################################
### 8.5.20) Total treatments required - Control
########################################
### Create a dichotomous variable
mods$total_treatments_required_cont_recode <- ifelse(is.na(mods$total_treatments_required_cont) == TRUE, 'FALSE', 'TRUE')
### Create an object for the row name, level_ref, and level_alternative names
mod_rep_treat_req_con_names <-c('reported_treatment_required_con','treatment', 'No', 'Yes')

### Create a vector/object of each subgroup's effect sizes, variance, pooled effect, Q and Q p_value
sg_rep_treat_req_con_No_es <-mods[mods$total_treatments_required_cont_recode== 'FALSE', "LRIRR"]
sg_rep_treat_req_con_No_v <-mods[mods$total_treatments_required_cont_recode== 'FALSE', "RIRR_var"]
sg_rep_treat_req_con_No_k <-length(sg_rep_treat_req_con_No_es)
sg_rep_treat_req_con_No_PLRIRR <-(sum(sg_rep_treat_req_con_No_es/sg_rep_treat_req_con_No_v))/(sum(1/sg_rep_treat_req_con_No_v))
sg_rep_treat_req_con_No_PLRIRR_v <- 1/sum(1/sg_rep_treat_req_con_No_v)
sg_rep_treat_req_con_No_Q <-sum((sg_rep_treat_req_con_No_es-sg_rep_treat_req_con_No_PLRIRR)^2/sg_rep_treat_req_con_No_v)
sg_rep_treat_req_con_No_Q_p <-pchisq(sg_rep_treat_req_con_No_Q, df=(length(sg_rep_treat_req_con_No_es)-1), lower.tail = FALSE)
sg_rep_treat_req_con_No_I2 <-((sg_rep_treat_req_con_No_Q-(length(sg_rep_treat_req_con_No_es)-1))/sg_rep_treat_req_con_No_Q)*100 

sg_rep_treat_req_con_Yes_es <-mods[mods$total_treatments_required_cont_recode== 'TRUE', "LRIRR"]
sg_rep_treat_req_con_Yes_v <-mods[mods$total_treatments_required_cont_recode== 'TRUE', "RIRR_var"]
sg_rep_treat_req_con_Yes_k <-length(sg_rep_treat_req_con_Yes_es)
sg_rep_treat_req_con_Yes_PLRIRR <-(sum(sg_rep_treat_req_con_Yes_es/sg_rep_treat_req_con_Yes_v))/(sum(1/sg_rep_treat_req_con_Yes_v))
sg_rep_treat_req_con_Yes_PLRIRR_v <- 1/sum(1/sg_rep_treat_req_con_Yes_v)
sg_rep_treat_req_con_Yes_Q <-sum((sg_rep_treat_req_con_Yes_es-sg_rep_treat_req_con_Yes_PLRIRR)^2/sg_rep_treat_req_con_Yes_v)
sg_rep_treat_req_con_Yes_Q_p <-pchisq(sg_rep_treat_req_con_Yes_Q, df=(length(sg_rep_treat_req_con_Yes_es)-1), lower.tail = FALSE)
sg_rep_treat_req_con_Yes_I2 <-((sg_rep_treat_req_con_Yes_Q-(length(sg_rep_treat_req_con_Yes_es)-1))/sg_rep_treat_req_con_Yes_Q)*100 

###Insert the moderator into the model
model.lm.mod.rep.treat.req.con<-lm(LRIRR ~ 1 + mods$other_treatment_treat_aftercare,weights = 1/v)

###Extract moderator Beta
mod_rep_treat_req_con_LRIRR<-summary(model.lm.mod.rep.treat.req.con)$coefficients[2, 1]

###Extract the Beta's standard error
mod_rep_treat_req_con_se<-summary(model.lm.mod.rep.treat.req.con)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_rep_treat_req_con_res_df<- model.lm.mod.rep.treat.req.con$df.residual

###Calculate Qres for the model
mod_rep_treat_req_con_Y <- as.matrix(mods$LRIRR)
mod_rep_treat_req_con_W <- diag(1/mods$RIRR_var)
mod_rep_treat_req_con_X <- model.matrix(model.lm.mod.rep.treat.req.con)
mod_rep_treat_req_con_P <- mod_rep_treat_req_con_W - mod_rep_treat_req_con_W %*% mod_rep_treat_req_con_X %*% solve(
  t(mod_rep_treat_req_con_X) %*% mod_rep_treat_req_con_W %*% mod_rep_treat_req_con_X) %*% t(
    mod_rep_treat_req_con_X) %*% mod_rep_treat_req_con_W 

mod_rep_treat_req_con_Qres <- max(0, c(crossprod(mod_rep_treat_req_con_Y, mod_rep_treat_req_con_P) %*% mod_rep_treat_req_con_Y))

###Bind the list of values to the dataframe 'moderator_output'
mod_rep_treat_req_con_list<-c(mod_rep_treat_req_con_names,sg_rep_treat_req_con_No_k, sg_rep_treat_req_con_No_PLRIRR,sg_rep_treat_req_con_No_PLRIRR_v, sg_rep_treat_req_con_No_Q, sg_rep_treat_req_con_No_Q_p, 
                             sg_rep_treat_req_con_No_I2,sg_rep_treat_req_con_Yes_k, sg_rep_treat_req_con_Yes_PLRIRR,sg_rep_treat_req_con_Yes_PLRIRR_v, sg_rep_treat_req_con_Yes_Q, 
                             sg_rep_treat_req_con_Yes_Q_p, sg_rep_treat_req_con_Yes_I2, mod_rep_treat_req_con_LRIRR, mod_rep_treat_req_con_se,
                             mod_rep_treat_req_con_res_df, mod_rep_treat_req_con_Qres)

names(mod_rep_treat_req_con_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                    "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                    "mod_Qres")

moderator_output<- bind_rows(moderator_output, mod_rep_treat_req_con_list)
########################################
### 8.5.21) Total treatments offered - Treatment
########################################
### Create an object for the row name, level_ref, and level_alternative names
tot_treat_off_tr_names<- c('total_treatments_offered_treat', 'treatment', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                           'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                           'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.tot.treat.off.tr<-lm(LRIRR ~ 1 + mods$total_treatments_offered_treat, weights = 1/v)

###Extract moderator Beta
mod_tot_treat_off_tr_B_LRIRR<-summary(model.lm.mod.tot.treat.off.tr)$coefficients[2, 1]

###Extract the Beta's standard error
mod_tot_treat_off_tr_B_se<-summary(model.lm.mod.tot.treat.off.tr)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_tot_treat_off_tr_res_df<- model.lm.mod.tot.treat.off.tr$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$total_treatments_offered_tr_is_na<- is.na(mods$total_treatments_offered_treat)
mod_tot_treat_off_tr_es <- mods[mods$total_treatments_offered_tr_is_na == 'FALSE', "LRIRR"]
mod_tot_treat_off_tr_v <- mods[mods$total_treatments_offered_tr_is_na == 'FALSE', "RIRR_var"]

mod_tot_treat_off_tr_Y <- as.matrix(mod_tot_treat_off_tr_es)
mod_tot_treat_off_tr_W <- diag(1/mod_tot_treat_off_tr_v)
mod_tot_treat_off_tr_X <- model.matrix(model.lm.mod.tot.treat.off.tr)
mod_tot_treat_off_tr_P <- mod_tot_treat_off_tr_W - mod_tot_treat_off_tr_W %*% mod_tot_treat_off_tr_X %*% solve(
  t(mod_tot_treat_off_tr_X) %*% mod_tot_treat_off_tr_W %*% mod_tot_treat_off_tr_X) %*% t(
    mod_tot_treat_off_tr_X) %*% mod_tot_treat_off_tr_W 

mod_tot_treat_off_tr_Qres <- max(0, c(crossprod(mod_tot_treat_off_tr_Y, mod_tot_treat_off_tr_P) %*% mod_tot_treat_off_tr_Y))

###Bind the list of values to the dataframe 'moderator_output'
tot_treat_off_tr_names_list<- c(tot_treat_off_tr_names, mod_tot_treat_off_tr_B_LRIRR, mod_tot_treat_off_tr_B_se, 
                                mod_tot_treat_off_tr_res_df, mod_tot_treat_off_tr_Qres)

names(tot_treat_off_tr_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                      "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                      "mod_Qres")
moderator_output<- bind_rows(moderator_output, tot_treat_off_tr_names_list)
########################################
### 8.5.22) Total treatments offered - control
########################################
### Create an object for the row name, level_ref, and level_alternative names
tot_treat_off_con_names<- c('total_treatments_offered_cont', 'treatment', 'Beta_weight', 'Not Applicable', 'Not Applicable', 'Not Applicable',  'Not Applicable', 'Not Applicable', 
                            'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable', 'Not Applicable',
                            'Not Applicable', 'Not Applicable', 'Not Applicable') 

###Insert the moderator into the model
model.lm.mod.tot.treat.off.con<-lm(LRIRR ~ 1 + mods$total_treatments_offered_cont, weights = 1/v)

###Extract moderator Beta
mod_tot_treat_off_con_B_LRIRR<-summary(model.lm.mod.tot.treat.off.con)$coefficients[2, 1]

###Extract the Beta's standard error
mod_tot_treat_off_con_B_se<-summary(model.lm.mod.tot.treat.off.con)$coefficients[2, 2]

###Create an object for the residual degrees of freedom. 
mod_tot_treat_off_con_res_df<- model.lm.mod.tot.treat.off.con$df.residual

###Calculate Qres for the model
###Separate values from NA
mods$total_treatments_offered_tr_is_na<- is.na(mods$total_treatments_offered_cont)
mod_tot_treat_off_con_es <- mods[mods$total_treatments_offered_tr_is_na == 'FALSE', "LRIRR"]
mod_tot_treat_off_con_v <- mods[mods$total_treatments_offered_tr_is_na == 'FALSE', "RIRR_var"]

mod_tot_treat_off_con_Y <- as.matrix(mod_tot_treat_off_con_es)
mod_tot_treat_off_con_W <- diag(1/mod_tot_treat_off_con_v)
mod_tot_treat_off_con_X <- model.matrix(model.lm.mod.tot.treat.off.con)
mod_tot_treat_off_con_P <- mod_tot_treat_off_con_W - mod_tot_treat_off_con_W %*% mod_tot_treat_off_con_X %*% solve(
  t(mod_tot_treat_off_con_X) %*% mod_tot_treat_off_con_W %*% mod_tot_treat_off_con_X) %*% t(
    mod_tot_treat_off_con_X) %*% mod_tot_treat_off_con_W 

mod_tot_treat_off_con_Qres <- max(0, c(crossprod(mod_tot_treat_off_con_Y, mod_tot_treat_off_con_P) %*% mod_tot_treat_off_con_Y))

###Bind the list of values to the dataframe 'moderator_output'
tot_treat_off_con_names_list<- c(tot_treat_off_con_names, mod_tot_treat_off_con_B_LRIRR, mod_tot_treat_off_con_B_se, 
                                 mod_tot_treat_off_con_res_df, mod_tot_treat_off_con_Qres)

names(tot_treat_off_con_names_list)<-c("name", "mod_category", "level_ref", "level_alt", "ref_sg_k", "ref_sg_PLRIRR", "ref_sg_PLRIRR_v", "ref_sg_Q", "ref_sg_Q_pv", "ref_sg_I2",
                                       "alt_sg_k", "alt_sg_PLRIRR", "alt_sg_PLRIRR_v","alt_sg_Q", "alt_sg_Q_pv", "alt_sg_I2", "Beta_LRIRR", "Beta_SE", "Res_DF", 
                                       "mod_Qres")
moderator_output<- bind_rows(moderator_output, tot_treat_off_con_names_list)

########################################
### 9) Z test for subgroup differences (categorical moderators) 
########################################
########################################
### 9.1) Multiplicative Variance Adjustment 
########################################
###Organise the rows into their moderator category and add an id to each row in moderator_output
moderator_output<- arrange(moderator_output, mod_category)
moderator_output$mod_id <- 1:nrow(moderator_output)
moderator_output<- moderator_output %>%
  select(mod_id, everything())
###Convert 'Not Applicables to NA
moderator_output <- na_if(moderator_output, 'Not Applicable')

###Higgins et al. (2003, p 558) note that negative values of I2 are equal to 0. 
###Convert any I2 values <0 to 0. 
moderator_output$ref_sg_I2<- ifelse(moderator_output$ref_sg_I2<0, 0, moderator_output$ref_sg_I2)
moderator_output$alt_sg_I2<- ifelse(moderator_output$alt_sg_I2<0, 0, moderator_output$alt_sg_I2)

###Convert relevant variables to numeric
moderator_output_num<-moderator_output%>%select(starts_with("ref"), starts_with("alt"),starts_with("beta"),
                                starts_with("mod"),-"mod_category", starts_with("res"))%>%
                                apply(2,function(x) as.numeric(as.character(x))) 

moderator_output_non_num<-moderator_output%>%select(-starts_with("ref"), -starts_with("alt"), -starts_with("beta"),
                                                -starts_with("mod"),"mod_category", -starts_with("res")) 


moderator_output_numbered<-cbind(moderator_output_non_num,moderator_output_num) #bind numeric and character variables
moderator_output_numbered<-moderator_output_numbered[,names(moderator_output)]

###Make the variance adjustment 
### Farrington & Welsh (2013) describe the following formula for the multiplicative variance adjustment
### Vm = Vf * Ve
### Where Vm = variance in MVA model
###       Vf = variance in the fixed effects model
###       Ve = extra variance where Ve = Q / df 

### MVA for the reference group 
moderator_output_numbered_ref_Ve <- moderator_output_numbered$ref_sg_Q/(moderator_output_numbered$ref_sg_k-1)
moderator_output_numbered$ref_Vm <- moderator_output_numbered$ref_sg_PLRIRR_v*moderator_output_numbered_ref_Ve

### MVA for the alt group
moderator_output_numbered_alt_Ve <- moderator_output_numbered$alt_sg_Q/(moderator_output_numbered$alt_sg_k-1)
moderator_output_numbered$alt_Vm <- moderator_output_numbered$alt_sg_PLRIRR_v*moderator_output_numbered_alt_Ve

### Calculate the test-statistic, standard error, test-statistic p-value, and 95% confidence interval for the ref and alt groups
moderator_output_numbered$ref_sg_MVA_z<-moderator_output_numbered$ref_sg_PLRIRR/sqrt(moderator_output_numbered$ref_Vm)
moderator_output_numbered$alt_sg_MVA_z<-moderator_output_numbered$alt_sg_PLRIRR/sqrt(moderator_output_numbered$alt_Vm)

moderator_output_numbered$ref_PLRIRR_se <-sqrt(moderator_output_numbered$ref_Vm) 
moderator_output_numbered$alt_PLRIRR_se <-sqrt(moderator_output_numbered$alt_Vm)

### P value
moderator_output_numbered$ref_z_p <- 2 * pnorm(abs(moderator_output_numbered$ref_sg_MVA_z), lower.tail = FALSE)
moderator_output_numbered$alt_z_p <- 2 * pnorm(abs(moderator_output_numbered$alt_sg_MVA_z), lower.tail = FALSE)

### 95% CI 
moderator_output_numbered$ref_z_CL <- exp(moderator_output_numbered$ref_sg_PLRIRR-1.96*sqrt(moderator_output_numbered$ref_Vm))
moderator_output_numbered$ref_z_CU <- exp(moderator_output_numbered$ref_sg_PLRIRR+1.96*sqrt(moderator_output_numbered$ref_Vm))

moderator_output_numbered$alt_z_CL <- exp(moderator_output_numbered$alt_sg_PLRIRR-1.96*sqrt(moderator_output_numbered$alt_Vm))
moderator_output_numbered$alt_z_CU <- exp(moderator_output_numbered$alt_sg_PLRIRR+1.96*sqrt(moderator_output_numbered$alt_Vm))

########################################
### 9.2) Z-test
########################################
### Perform the 'Z-test method' for differences between groups described by Borenstein et al (2009, p156):
### Diff = Mb - Ma
### ZDiff = Diff/ SEdiff
### SEdiff = sqrt(Vma + Vmb) 
###                          Where Mb and Ma are the two pooled effects, 
###                          and Vma and Vmb are the variance of the two pooled effects 

### Calculate Diff, SEdiff, and ZDiff
moderator_output_numbered$Diff <- (moderator_output_numbered$alt_sg_PLRIRR - moderator_output_numbered$ref_sg_PLRIRR)
moderator_output_numbered$SEdiff <- sqrt(moderator_output_numbered$ref_Vm + moderator_output_numbered$alt_Vm)
moderator_output_numbered$Zdiff <- moderator_output_numbered$Diff/moderator_output_numbered$SEdiff

### A a two tailed p-value can then be calculated on the ZDiff using a normal distribution (Borenstein et al., 2009, p. 156)
moderator_output_numbered$Zdiff_p <- 2 * pnorm(abs(moderator_output_numbered$Zdiff), lower.tail = FALSE)

########################################
### 10) Benjamini-Hochberg  adjustments  
########################################
########################################
### 10.1) Categorical moderators
########################################
###First seperate the categorical moderators from the continuous 
moderator_output_cat <- moderator_output_numbered%>%filter(!level_ref=="Beta_weight")

###Check that each covariate has at least k = 10
moderator_k_total <- aggregate(moderator_output_cat$alt_sg_k, by=list(Category=moderator_output_cat$name), FUN=sum)
moderator_k_total$ref <- moderator_output_cat$ref_sg_k[match(moderator_k_total$Category,moderator_output_cat$name)]
moderator_k_total$tot_k <- moderator_k_total$x + moderator_k_total$ref
moderator_k_total # each covariate has more than k = 10 studies, and each subgroup has more than k = 2. 

###Seperate a category from moderator analysis
mod_category_methodological <-moderator_output_cat%>%filter(mod_category=='methodological')
###Create a new variable that for BH adjusted p-values 
mod_category_methodological$Zdiff_p_adjusted <- p.adjust(mod_category_methodological$Zdiff_p, 
                                  method = "BH")
###Repeat for the remaining moderator categories 
###bias_&_quality
mod_category_quality_and_bias <-moderator_output_cat%>%filter(mod_category=='Quality_and_Bias')
mod_category_quality_and_bias$Zdiff_p_adjusted <- p.adjust(mod_category_quality_and_bias$Zdiff_p, 
                                                               method = "BH")
###sample_characteristics 
mod_category_sample_characteristics <-moderator_output_cat%>%filter(mod_category=='sample_characteristics')
mod_category_sample_characteristics$Zdiff_p_adjusted <-p.adjust(mod_category_sample_characteristics$Zdiff_p, 
                                                                     method = "BH")
###programmatic
mod_category_programmatic<-moderator_output_cat%>%filter(mod_category=='programmatic')
mod_category_programmatic$Zdiff_p_adjusted <-p.adjust(mod_category_programmatic$Zdiff_p, 
                                                                     method = "BH")
###treatment
mod_category_treatment<-moderator_output_cat%>%filter(mod_category=='treatment')
mod_category_treatment$Zdiff_p_adjusted <-p.adjust(mod_category_treatment$Zdiff_p, 
                                                          method = "BH")
###bind the dataframes together 
moderator_output_cat_adjusted <-rbind(mod_category_methodological,
                                  mod_category_quality_and_bias,
                                  mod_category_sample_characteristics,
                                  mod_category_programmatic,
                                  mod_category_treatment)
###Re-order by mod_id
moderator_output_cat_adjusted<- arrange(moderator_output_cat_adjusted, mod_id)
###check that the values from moderator_output_cat and moderator_output_cat_adjusted are the same 
moderator_output_cat$mod_id==moderator_output_cat_adjusted$mod_id
moderator_output_cat$Zdiff_p==moderator_output_cat_adjusted$Zdiff_p
###add a column for rounded p-values 
moderator_output_cat_adjusted$Zdiff_p_adjusted_rounded <-round(moderator_output_cat_adjusted$Zdiff_p_adjusted, digits = 3)

########################################
### 10.2) Continuous moderators
########################################
### First separate the categorical moderators from the continuous 
moderator_output_cont <- moderator_output_numbered%>%filter(level_ref=="Beta_weight")

### Add a column with each moderator's k
moderator_output_cont$beta_k <- moderator_output_cont$Res_DF + 2

### Deeks et al., (2008, p. 284) note that there should be 10 studies for each covariate.
### Remove any moderator that has fewer than k = 10 studies per outcome. 
moderator_output_cont <- moderator_output_cont[c(moderator_output_cont$beta_k >= 10), ]

### Create a new variable for the Z statistic following Borenstein et al.'s (2009, p 191) formula
moderator_output_cont$beta_z <- moderator_output_cont$Beta_LRIRR/ moderator_output_cont$Beta_SE

### Create variables for the Z statistic's p value and 95% confidence interval
moderator_output_cont$beta_z_p <- 2 * pnorm(abs(moderator_output_cont$beta_z), lower.tail = FALSE)

moderator_output_cont$beta_z_CL <- exp(moderator_output_cont$Beta_LRIRR-1.96*moderator_output_cont$Beta_SE)
moderator_output_cont$beta_z_CU <- exp(moderator_output_cont$Beta_LRIRR+1.96*moderator_output_cont$Beta_SE)
### Create a new variable that for BH adjusted p-values 
moderator_output_cont$beta_p_adjusted <- p.adjust(moderator_output_cont$beta_z_p, 
                                                         method = "BH")

### Select relevant variables
moderator_output_cont_adjusted <- select(moderator_output_cont, -level_ref,-level_alt, -starts_with("ref"), -starts_with("alt"),-ends_with("Vm"), -beta_z, -ends_with("CU"),
                                      -ends_with("CL"),beta_k, beta_z_CL, beta_z_CU, -Zdiff,-Zdiff_p, -SEdiff, -Diff, -ref_sg_MVA_z, -alt_sg_MVA_z,
                                      -ref_z_p, -alt_z_p)
########################################
### 11) Pooled effects table 
########################################
###Convert relevant columns to numeric
pooled_effect_col_names<- c("Relative Incident Rate Ratio", "Conlow", "Conup", "Test statistic", "p-value",
                             "credlow", "credup", "Q statistic", "Q_pv", "I2")
pooled_effect_output[pooled_effect_col_names] <- sapply(pooled_effect_output[pooled_effect_col_names], as.numeric) #ignore following warning message
sapply(pooled_effect_output, class)
###Convert p-values <.001 to <.001
pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {
  
  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }
  
  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit)
      if (html)
        return(sprintf('&lt; %s', format(sig.limit))) else
          return(sprintf('< %s', format(sig.limit)))
    if (x > .1)
      return(roundr(x, digits = 2)) else
        return(roundr(x, digits = digits))
  }, sig.limit = sig.limit)
}

pooled_effect_output <-pooled_effect_output %>% mutate_at(vars(`p-value`, Q_pv), 
                                                          list(~pvalr(., digits = 3)))
### Round columns to two decimal places 
pooled_effect_output <- pooled_effect_output %>% mutate_at(vars(-Model, -`p-value`, -Q_pv), list(~round(., 2)))

###Add trailing zeros to the CI numbers 
pooled_effect_output <-pooled_effect_output %>% mutate_at(vars(Conlow, Conup, `Test statistic`, credlow, credup, I2), 
                                                          list(~sprintf("%.2f",. )))
###Join the interval variables, place them inside brackests, and place a comma between them
pooled_effect_output <-unite(pooled_effect_output , 'Confidence_Interval', c(Conlow, Conup), remove=TRUE)
pooled_effect_output <-unite(pooled_effect_output , 'Credibility_Interval', c(credlow, credup), remove=TRUE)

pooled_effect_output$`Confidence_Interval` <-  paste0("(", (pooled_effect_output$`Confidence_Interval`), ")")
pooled_effect_output$`Credibility_Interval` <-  paste0("(", (pooled_effect_output$`Credibility_Interval`), ")")

pooled_effect_output$`Confidence_Interval` <-gsub("_", ", ", pooled_effect_output$`Confidence_Interval`)
pooled_effect_output$`Credibility_Interval` <-gsub("_", ", ", pooled_effect_output$`Credibility_Interval`)
pooled_effect_output$`Credibility_Interval` <- ifelse(pooled_effect_output$`Credibility_Interval` == "(NA, NA)", ' - ', pooled_effect_output$`Credibility_Interval`)

###Add ***to reflect the significance level of p-values
pooled_effect_output$`p-value`<- ifelse(pooled_effect_output$`p-value`== "< 0.001", paste0(pooled_effect_output$`p-value`,"***"),
                                              pooled_effect_output$`p-value`)

pooled_effect_output$Q_pv<- ifelse(pooled_effect_output$Q_pv== "< 0.001", paste0(pooled_effect_output$Q_pv,"***"),
                                        pooled_effect_output$Q_pv)

# write.csv(pooled_effect_output, file = "pooled_effect_output_table.csv")

########################################
### 12) Categorical moderator table 
########################################
########################################
### 12.1) Format table
########################################
###First separate the reference and alternative subgroups into different rows
moderator_output_ref <- select(moderator_output_cat_adjusted, -ends_with("Vm"), -starts_with("Beta"), 
                               -starts_with("alt"), -ends_with("alt"), -ends_with("PLRIRR_v"), -Res_DF, -mod_Qres, -Zdiff, 
                               -SEdiff, -ref_sg_MVA_z, -Zdiff_p_adjusted)
###Replace cells that will later be formatted blank with NA
moderator_output_ref$Diff<- NA
moderator_output_ref$Zdiff_p<- NA
moderator_output_ref$Zdiff_p_adjusted_rounded<- NA

###create a variable stating whether the row is ref/alt
moderator_output_ref$group <- "ref"

###Remove ref from variable names 
names(moderator_output_ref) <- sub("ref_", "", names(moderator_output_ref))
names(moderator_output_ref) <- sub("_ref", "", names(moderator_output_ref))

###Additionally convert relevant p-values to <.001
moderator_output_ref <-moderator_output_ref %>% mutate_at(vars(sg_Q_pv, z_p), 
                                                          list(~pvalr(., digits = 3)))
###Repeat with the alt subgroups 
moderator_output_alt <- select(moderator_output_cat_adjusted, -ends_with("Vm"), -starts_with("Beta"), 
                               -starts_with("ref"), -ends_with("ref"), -ends_with("PLRIRR_v"), -Res_DF, -mod_Qres, -Zdiff, 
                               -SEdiff, -alt_sg_MVA_z, -Zdiff_p_adjusted)
moderator_output_alt$group <- "alt"
names(moderator_output_alt) <- sub("alt_", "", names(moderator_output_alt))
names(moderator_output_alt) <- sub("_alt", "", names(moderator_output_alt))

moderator_output_alt <-moderator_output_alt %>% mutate_at(vars(sg_Q_pv, z_p, Zdiff_p, 
                                                                           Zdiff_p_adjusted_rounded), 
                                                                                list(~pvalr(., digits = 3)))

###bind the two dataframes together, group them by mod id - need to exclude rows where there is more than one alt
moderator_output_cat_final <- arrange(rbind(moderator_output_ref,moderator_output_alt), mod_id)
moderator_output_cat_final <- arrange(moderator_output_cat_final, mod_id)

###Remove excess rows where there are multiple alternative subrgroups
moderator_output_cat_final <-distinct(moderator_output_cat_final, name, level, .keep_all= TRUE)

###Exponentiate relevant columns and change variable names accordingly
moderator_output_cat_final <-moderator_output_cat_final %>% mutate_at(vars(sg_PLRIRR, Diff),
                                                                                 list(~exp(.)))
moderator_output_cat_final<- moderator_output_cat_final %>% rename(RES = sg_PLRIRR, Subgroups_difference = Diff)


###Round relevant columns 
moderator_output_cat_final <-moderator_output_cat_final %>% mutate_at(vars(RES, sg_Q, sg_I2, 
                                                                           z_CL, z_CU, Subgroups_difference, 
                                                                           PLRIRR_se), list(~round(., 2)))
########################################
### 12.2) Add text, move and merge columns, and export 
########################################
###Re-order 
moderator_output_cat_final<- moderator_output_cat_final[ , c('mod_category', 'name', 'level', 'sg_k','RES','z_p','PLRIRR_se','z_CL','z_CU','sg_Q','sg_Q_pv','sg_I2',
                                                            'Subgroups_difference', 'Zdiff_p','Zdiff_p_adjusted_rounded','mod_id', 'group')]

###Paste in text
moderator_output_cat_final$level<- ifelse((moderator_output_cat_final$group== 'ref') == TRUE, paste0(moderator_output_cat_final$level," (reference, k = "), 
                                          moderator_output_cat_final$level)
moderator_output_cat_final$level<- ifelse((moderator_output_cat_final$group== 'alt') == TRUE, paste0(moderator_output_cat_final$level," (k = "), 
                                          moderator_output_cat_final$level)
moderator_output_cat_final$sg_k<- paste0(moderator_output_cat_final$sg_k,")")

###Join level and sg_k variables 
moderator_output_cat_final <-unite(moderator_output_cat_final , 'Subgroup', c(level, sg_k), remove=TRUE)
moderator_output_cat_final$Subgroup <- str_replace(moderator_output_cat_final$Subgroup, "_", "")

###Add trailing zeros to the CI numbers 
moderator_output_cat_final <-moderator_output_cat_final %>% mutate_at(vars(z_CL, z_CU), 
                                                                      list(~sprintf("%.2f",. )))

###Add significant asterisks 
moderator_output_cat_final$sg_Q<- ifelse(moderator_output_cat_final$sg_Q_pv > 0.01 & moderator_output_cat_final$sg_Q_pv<0.05, paste0(moderator_output_cat_final$sg_Q,"*"),
                                            ifelse(moderator_output_cat_final$sg_Q_pv>0.001 & moderator_output_cat_final$sg_Q_pv<0.01, paste0(moderator_output_cat_final$sg_Q,"**"),
                                                   ifelse(moderator_output_cat_final$sg_Q_pv< 0.001, paste0(moderator_output_cat_final$sg_Q,"***"), 
                                                          moderator_output_cat_final$sg_Q)))

moderator_output_cat_final$Zdiff_p<- ifelse(moderator_output_cat_final$Zdiff_p_adjusted_rounded > 0.01 & moderator_output_cat_final$Zdiff_p_adjusted_rounded<0.05, paste0(moderator_output_cat_final$Zdiff_p,"*"),
                                            ifelse(moderator_output_cat_final$Zdiff_p_adjusted_rounded>0.001 & moderator_output_cat_final$Zdiff_p_adjusted_rounded<0.01, paste0(moderator_output_cat_final$Zdiff_p,"**"),
                                                   ifelse(moderator_output_cat_final$Zdiff_p_adjusted_rounded< 0.001, paste0(moderator_output_cat_final$Zdiff_p,"***"), 
                                                          moderator_output_cat_final$Zdiff_p)))

###Join the interval variables, place them inside brackets, and place a comma between them
moderator_output_cat_final <-unite(moderator_output_cat_final , 'z_CI', c(z_CL, z_CU), remove=TRUE)

moderator_output_cat_final$z_CI <-paste0("(", (moderator_output_cat_final$z_CI), ")")
moderator_output_cat_final$z_CI <-gsub("_", ", ", moderator_output_cat_final$z_CI)

###replace unwanted strings 
moderator_output_cat_final$name <- stri_replace_last_regex(moderator_output_cat_final$name, '_', '-')
moderator_output_cat_final <- moderator_output_cat_final %>%
                                        separate(name, c("name_temp", "group_name_replace"), "-") # warning message can be ignored

moderator_output_cat_final$group_name_replace <- gsub("treat", "(T)", moderator_output_cat_final$group_name_replace)
moderator_output_cat_final$group_name_replace <- gsub("cont", "(C)", moderator_output_cat_final$group_name_replace)
moderator_output_cat_final <-unite(moderator_output_cat_final , 'Moderator', c(name_temp, group_name_replace), remove=TRUE)
moderator_output_cat_final$Moderator <- gsub("_NA", "", moderator_output_cat_final$Moderator)
moderator_output_cat_final$Moderator <- gsub("_", " ", moderator_output_cat_final$Moderator)

###Capitalise the first letter of $Moderator
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

moderator_output_cat_final <- moderator_output_cat_final %>% mutate_at(vars(Moderator, Subgroup), 
                                                                       list(~firstup(.)))

###Replace duplicate 'moderator' values 
moderator_output_cat_final$Moderator[duplicated(moderator_output_cat_final$Moderator)] <- NA

###Select the final variables to export and export
moderator_output_cat_final<- moderator_output_cat_final[, c("mod_category", "Moderator", "Subgroup", "RES",                 
                                                                        "PLRIRR_se", "z_CI", "sg_Q", "sg_I2", 
                                                                        "Subgroups_difference", "Zdiff_p")]

# write.csv(moderator_output_cat_final, file = "moderator_table_categorical.csv", na = "")

########################################
### 13) Continuous moderator table
########################################
########################################
### 13.1) Amend figures
########################################
### Calculate Qres' p-value and the model I2 
moderator_output_cont_adjusted$Qres_pv <-pvalr(pchisq(moderator_output_cont_adjusted$mod_Qres, 
                                                df=(moderator_output_cont_adjusted$Res_DF), lower.tail = FALSE), 
                                                digits = 3)
moderator_output_cont_adjusted$I2 <-((moderator_output_cont_adjusted$mod_Qres-(moderator_output_cont_adjusted$Res_DF))/moderator_output_cont_adjusted$mod_Qres)*100  

### Reorder
moderator_output_cont_final<- moderator_output_cont_adjusted[ , c("mod_category", "name", "beta_k", "Beta_LRIRR","beta_z_CL", "beta_z_CU", "beta_p_adjusted", "Beta_SE", 
                                                                  "mod_Qres", "Qres_pv", "Res_DF", "I2", "mod_id")]
###Exponentiate relevant columns and change variable names
moderator_output_cont_final$Beta_LRIRR <-exp(moderator_output_cont_final$Beta_LRIRR)
moderator_output_cont_final <-moderator_output_cont_final %>% rename(b = Beta_LRIRR, SE = Beta_SE, Qres = mod_Qres, 
                                                                     Residual_df = Res_DF, k = beta_k)
###Round relevant columns and add trailing 0s 
no_one <- function(x){ifelse(x >= .995 & x < 1, round(x, 3), (sprintf("%.2f", round(x, 2))))
  }

moderator_output_cont_final <-moderator_output_cont_final %>% mutate_at(vars(b, beta_z_CL, beta_z_CU, SE, Qres, I2),  list(~no_one(.)))
                                                              
###Add asterisks to p-values
moderator_output_cont_final$Qres<- ifelse(moderator_output_cont_final$Qres_pv > 0.01 & moderator_output_cont_final$Qres_pv<0.05, paste0(moderator_output_cont_final$Qres,"*"),
                                        ifelse(moderator_output_cont_final$Qres_pv>0.001 & moderator_output_cont_final$Qres_pv<0.01, paste0(moderator_output_cont_final$Qres,"**"),
                                               ifelse(moderator_output_cont_final$Qres_pv< 0.001, paste0(moderator_output_cont_final$Qres,"***"), 
                                                      moderator_output_cont_final$Qres)))

moderator_output_cont_final$b<- ifelse(moderator_output_cont_final$beta_p_adjusted > 0.01 & moderator_output_cont_final$beta_p_adjusted<0.05, paste0(moderator_output_cont_final$b,"*"),
                                              ifelse(moderator_output_cont_final$beta_p_adjusted>0.001 & moderator_output_cont_final$beta_p_adjusted<0.01, paste0(moderator_output_cont_final$b,"**"),
                                                     ifelse(moderator_output_cont_final$beta_p_adjusted< 0.001, paste0(moderator_output_cont_final$b,"***"), 
                                                            moderator_output_cont_final$b)))
########################################
### 13.2) Edit text and export
########################################
###Join the interval variables, place them inside brackests, and place a comma between them
moderator_output_cont_final <-unite(moderator_output_cont_final , 'CI', c(beta_z_CL, beta_z_CU), remove=TRUE)

moderator_output_cont_final$CI <-paste0("(", (moderator_output_cont_final$CI), ")")

moderator_output_cont_final$CI <-gsub("_", ", ", moderator_output_cont_final$CI)
###replace unwanted strings 
moderator_output_cont_final$name<- stri_replace_last_regex(moderator_output_cont_final$name, '_', '/')
moderator_output_cont_final <- moderator_output_cont_final %>%
                                        separate(name, c("name_temp", "group_name_replace"), "/")

moderator_output_cont_final$group_name_replace <- gsub("treat", "(T)", moderator_output_cont_final$group_name_replace)
moderator_output_cont_final$group_name_replace <- gsub("cont", "(C)", moderator_output_cont_final$group_name_replace)
moderator_output_cont_final <-unite(moderator_output_cont_final , 'Moderator', c(name_temp, group_name_replace), remove=TRUE)
moderator_output_cont_final$Moderator <- gsub("_", " ", moderator_output_cont_final$Moderator)

###capitalise first letter of $Moderator
moderator_output_cont_final$Moderator <- firstup(moderator_output_cont_final$Moderator)

###Drop redundant columns
moderator_output_cont_final <- moderator_output_cont_final[ , c('mod_category', 'Moderator','k', 'b', 'CI','SE',
                                                                'Qres','I2')]
                                                                
###Export table 
# write.csv(moderator_output_cont_final, file = "moderator_table_continuous.csv", na = "")

########################################
### 14) Health and well-being outcomes
########################################
########################################
### 14.1) Separate target rows 
########################################
health_out <- subset(MA_EX_PT_final,(Outcome=='Health'))

### Separate rows to preferred measures
health_meas <- subset(health_out, subset = Measure %in% c("crisis services received", 
                                                              "mean CMH substance abuse services",
                                                              "mean CFS substance abuse services", 
                                                              "mean crisis episodes",
                                                              "mean high intensity services received",
                                                              "mean intensive treatment episodes", 
                                                              "mean low intensity services received",
                                                              "mean medium intensity services received",
                                                              "mean therapeutic treatment episodes",
                                                              "non-crisis services received"))

health_pref_ITT<- subset(health_meas, subset = ITT %in% "Graduates/Failures")
health_unpref_ITT <- subset(health_meas, !(health_meas$ITT== "Graduates/Failures")) 
health_unpref_ITT <- health_unpref_ITT[!(health_unpref_ITT$studyid%in%health_pref_ITT$studyid),] 
health_bound<-bind_rows(health_pref_ITT, health_unpref_ITT) 

sum(unique(health_meas$studyid)>0)==sum(unique(health_bound$studyid)>0)

########################################
### 14.2) Calculate effect sizes and statistical significance 
########################################
### Note: we will calculate the RIRR following the same process outlined in Section '3) Calculate the Relative Incident Rate Ratio'.
### Calculate each RIRR
health_bound$RIRR <- ifelse(is.na(health_bound$bl_treat_m & health_bound$bl_cont_m & health_bound$pt_treat_m & health_bound$pt_cont_m) == TRUE,
                          ((health_bound$bl_treat_Dichotomous_._Yes*health_bound$pt_cont_Dichotomous_.))/
                          ((health_bound$pt_treat_Dichotomous_.*health_bound$bl_cont_Dichotomous_._Yes)),
                          ((health_bound$bl_treat_m*health_bound$pt_cont_m)/
                          (health_bound$pt_treat_m*health_bound$bl_cont_m)))

### Calculate each RIRR's variance 
health_bound$RIRR_var <- ifelse(is.na(health_bound$bl_treat_m & health_bound$bl_cont_m & health_bound$pt_treat_m & health_bound$pt_cont_m) == TRUE,
                               ((1 / health_bound$bl_treat_Dichotomous_._Yes) + 
                               (1 / health_bound$pt_treat_Dichotomous_.) +
                               (1 / health_bound$bl_cont_Dichotomous_._Yes) + 
                               (1 / health_bound$pt_cont_Dichotomous_.)), 
                               ((1 / health_bound$bl_in_rounded_treat) + 
                               (1 / health_bound$pt_in_rounded_treat) +
                               (1 / health_bound$bl_in_rounded_cont) + 
                               (1 / health_bound$pt_in_rounded_cont))) 

table(health_bound$Measure, round(health_bound$RIRR, digits = 3))

### Calculate a probability value 
### Z score: RES
health_bound$RIRR_z <-log(health_bound$RIRR)/sqrt(health_bound$RIRR_var)
### Z score: p value
health_bound$RIRR_p <- 1-pnorm(health_bound$RIRR_z)  
### Z score confidence interval 
health_bound$RIRR_z_cl<- exp(log(health_bound$RIRR)-1.96*sqrt(health_bound$RIRR_var))
health_bound$RIRR_z_cu<- exp(log(health_bound$RIRR)+1.96*sqrt(health_bound$RIRR_var))
### Round each to 2 decimal places (except p: round this to 3)
table(health_bound$Measure, pvalr(health_bound$RIRR_p, digits = 3))