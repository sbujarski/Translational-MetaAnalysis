#Translational Meta-Analysis R Code
#Goals of this study are to 
# 1) Collect data on all published trials using either:
#   A) Randomized Controlled Trial of a medication for an Alcohol Use Disorder (AUD)
#   B) Laboratory study testing the effects of a medication on responses to alcohol administration
# 2) Test the concordance between these different outcomes across the ranges of medications tested

library(SpPack) #Custom Function Suite
library(ggplot2) #Graphing
library(dplyr) #Data Wrangling
library(MAd) #agg function Implements Borenstein et al. (2009) approach for combining dependent ES and ES variances. Use Hedge's G
library(metafor) #meta-regression pachage
library(readxl) #package to import excel files directly
library(compute.es)


#CUSTOM FUNCTIONS------

#gg.funnel
#produces a funnel plot with effect sizes and variances
#medication labeling optional with lab (default = NULL)
gg.funnel <- function(es, es.var, mean.effect, se.effect, title, x.lab, y.lab, lab=NULL, labsTitle="gg.funnel")
{
  if(!is.null(lab))
  {
    rawdata<-data.frame(es=es, es.se=sqrt(es.var), lab=lab)
    
    se.max = ceiling(max(rawdata$es.se)*10)/10
    lwr95.limit <- mean.effect-1.96*se.max
    upr95.limit <- mean.effect+1.96*se.max
    CI95line=data.frame(x=c(lwr95.limit,mean.effect,upr95.limit),
                        y=c(se.max,0,se.max))
    CI95shade=data.frame(x=c(lwr95.limit, lwr95.limit, upr95.limit, upr95.limit, mean.effect, lwr95.limit),
                         y=c(se.max,      0,           0,           se.max,      0,           se.max))
    
    ggplot(data=rawdata, aes(x=es, y=es.se))+
      geom_line(data=CI95line, aes(x=x, y=y), linetype="dashed", size=1)+
      geom_polygon(data=CI95shade, aes(x=x, y=y), alpha=.1, fill="black")+
      geom_vline(xintercept=mean.effect, size=1)+
      geom_vline(xintercept=0, linetype="dotted", size=1)+
      geom_point(size=3, aes(colour=lab))+
      geom_text(size=3, hjust = 0, nudge_x = 0.03, check_overlap = TRUE, aes(label=abbreviate(lab, minlength = 3)))+
      ggtitle(title)+
      xlab(x.lab)+
      ylab(y.lab)+
      SpTheme()+
      scale_y_reverse()+
      guides(colour=guide_legend(title=labsTitle))
  }
  else
  {
    rawdata<-data.frame(es=es, es.se=sqrt(es.var))
    
    se.max = ceiling(max(rawdata$es.se)*10)/10
    lwr95.limit <- mean.effect-1.96*se.max
    upr95.limit <- mean.effect+1.96*se.max
    CI95line=data.frame(x=c(lwr95.limit,mean.effect,upr95.limit),
                        y=c(se.max,0,se.max))
    CI95shade=data.frame(x=c(lwr95.limit, lwr95.limit, upr95.limit, upr95.limit, mean.effect, lwr95.limit),
                         y=c(se.max,      0,           0,           se.max,      0,           se.max))
    
    ggplot(data=rawdata, aes(x=es, y=es.se))+
      geom_line(data=CI95line, aes(x=x, y=y), linetype="dashed", size=1)+
      geom_polygon(data=CI95shade, aes(x=x, y=y), alpha=.1, fill="black")+
      geom_vline(xintercept=mean.effect, size=1)+
      geom_vline(xintercept=0, linetype="dotted", size=1)+
      geom_point(size=3)+
      ggtitle(title)+
      xlab(x.lab)+
      ylab(y.lab)+
      SpTheme()+
      scale_y_reverse()
  }
}

#rma.recenterMed.Lab
#Function takes a pre-processed dataset with Med, es, var and covariate values for Laboratory data
#runs random effects meta-analysis with med, MaxDose.C (log-base2 modal centered medication dose), MaxAlcDose.C (0.06 centered), 
#                                       and Pop (population from light drinking, heavy drinking, or AUD; centered at AUD)
#prints meta analysis results with first med as reference group
#prints overall average effect size and se from a separate meta-analysis without med predictor for use in gg.funnel plot
#systematically recenters med factor and run rma to get estimate of es (intercept when centered) for each medication
#prints metaES, se, pval, ci.lb, and ci.ub for each medication based on meta-regression model
#returns uncentered rma result, recentered effect size estimates dataframe, overall mean, and overall SEM
rma.recenterMed.Lab <- function(data, abr=NULL)
{
  nMeds <- length(levels(data$Med))
  #print initial results 
  contrasts(data$Med) <- contr.treatment(nMeds, base=1)
  contrasts(data$Pop) <- contr.treatment(3, base=3)
  cat("Random Effects Meta-Analysis result", "\n")
  rma.uncent <- rma(yi=es, vi=var, mods = ~ Med  + MaxDose.C + MaxAlcDose.C + Pop, #meta regression outcome
                    data = data, method="REML")
  print(rma.uncent)
  
  #Get results averaging accross meds for funnel plot
  rma.noMed <- rma(yi=es, vi=var, mods = ~ MaxDose.C + MaxAlcDose.C + Pop, data = data, method="REML")
  cat(c("Overall Mean: ", as.double(rma.noMed$b["intrcpt",1])), "\n")
  cat(c("Overall SEM:  ", as.double(rma.noMed$se[1])), "\n", "\n")
  
  #Set up data frame
  ES.est <- data.frame(Med=levels(data$Med), metaES = rep(NA,nMeds), metaES.se = rep(NA,nMeds), metaES.pval = rep(NA,nMeds),
                       metaES.ci.lb = rep(NA,nMeds), metaES.ci.ub = rep(NA,nMeds))
  
  
  #Start contrasts
  for(i in 1:nMeds)
  {
    contrasts(data$Med) <- contr.treatment(nMeds, base=i)
    #print(contrasts(data$Med))
    
    #Run Meta-analysis
    rma.result <- rma(yi=es, vi=var, mods = ~ Med  + MaxDose.C + MaxAlcDose.C + Pop, #meta regression outcome
                      data = data, method="REML")
    #print(rma.result)
    ES.est$metaES[i] <- rma.result$b["intrcpt",1]
    ES.est$metaES.se[i] <- rma.result$se[1]
    ES.est$metaES.pval[i] <- rma.result$pval[1]
    ES.est$metaES.ci.lb[i] <- rma.result$ci.lb[1]
    ES.est$metaES.ci.ub[i] <- rma.result$ci.ub[1]
  }
  
  if(is.null(abr))
  {
    cat("Recentered Effect Size Estimates", "\n")
    print(ES.est)
    return(list(rma.uncent=rma.uncent, ES.est=ES.est, 
                ES.mean=as.double(rma.noMed$b["intrcpt",1]), ES.SEM=as.double(rma.noMed$se[1])))
  }
  else
  {
    names(ES.est)[names(ES.est)=="metaES"] <- paste(abr,"metaES",sep="")
    names(ES.est)[names(ES.est)=="metaES.se"] <- paste(abr,"metaES.se",sep="")
    names(ES.est)[names(ES.est)=="metaES.pval"] <- paste(abr,"metaES.pval",sep="")
    names(ES.est)[names(ES.est)=="metaES.ci.lb"] <- paste(abr,"metaES.ci.lb ",sep="")
    names(ES.est)[names(ES.est)=="metaES.ci.ub"] <- paste(abr,"metaES.ci.ub",sep="")
    
    cat("Recentered Effect Size Estimates", "\n")
    print(ES.est)
    return(list(rma.uncent=rma.uncent, ES.est=ES.est, 
                ES.mean=as.double(rma.noMed$b["intrcpt",1]), ES.SEM=as.double(rma.noMed$se[1])))
  }
}

#rma.recenterMed.RCT
#Function takes a pre-processed dataset with Med, es, var and covariate values for RCT data
#runs random effects meta-analysis with med, MaxDose.C (log-base2 modal centered medication dose), 
#                                       and TrxDur (treatment duration in weeks centered at 12)
#prints meta analysis results with first med as reference group
#prints overall average effect size and se from a separate meta-analysis without med predictor for use in gg.funnel plot
#systematically recenters med factor and run rma to get estimate of es (intercept when centered) for each medication
#prints metaES, se, pval, ci.lb, and ci.ub for each medication based on meta-regression model
#returns uncentered rma result, recentered effect size estimates dataframe, overall mean, and overall SEM
rma.recenterMed.RCT <- function(data, abr=NULL)
{
  nMeds <- length(levels(data$Med))
  #print initial results 
  contrasts(data$Med) <- contr.treatment(nMeds, base=1)
  
  #run first rma
  cat("Random Effects Meta-Analysis result", "\n")
  rma.uncent <- rma(yi=es, vi=var, mods = ~ Med  + MaxDose.C + TrxDur.C, #meta regression outcome
                    data = data, method="REML")
  print(rma.uncent)
  
  #Get results averaging accross meds for funnel plot
  rma.noMed <- rma(yi=es, vi=var, mods = ~ MaxDose.C + TrxDur.C, data = data, method="REML")
  cat(c("Overall Mean: ", as.double(rma.noMed$b["intrcpt",1])), "\n")
  cat(c("Overall SEM:  ", as.double(rma.noMed$se[1])), "\n", "\n")
  
  #Set up data frame
  ES.est <- data.frame(Med=levels(data$Med), metaES = rep(NA,nMeds), metaES.se = rep(NA,nMeds), metaES.pval = rep(NA,nMeds),
                       metaES.ci.lb = rep(NA,nMeds), metaES.ci.ub = rep(NA,nMeds))
  
  
  #Start contrasts
  for(i in 1:nMeds)
  {
    contrasts(data$Med) <- contr.treatment(nMeds, base=i)
    #print(contrasts(data$Med))
    
    #Run Meta-analysis
    rma.result <- rma(yi=es, vi=var, mods = ~ Med  + MaxDose.C + TrxDur.C, #meta regression outcome
                      data = data, method="REML")
    #print(rma.result)
    ES.est$metaES[i] <- rma.result$b["intrcpt",1]
    ES.est$metaES.se[i] <- rma.result$se[1]
    ES.est$metaES.pval[i] <- rma.result$pval[1]
    ES.est$metaES.ci.lb[i] <- rma.result$ci.lb[1]
    ES.est$metaES.ci.ub[i] <- rma.result$ci.ub[1]
  }
  
  if(is.null(abr))
  {
    cat("Recentered Effect Size Estimates", "\n")
    print(ES.est)
    return(list(rma.uncent=rma.uncent, ES.est=ES.est, 
                ES.mean=as.double(rma.noMed$b["intrcpt",1]), ES.SEM=as.double(rma.noMed$se[1])))
  }
  else
  {
    names(ES.est)[names(ES.est)=="metaES"] <- paste(abr,"metaES",sep="")
    names(ES.est)[names(ES.est)=="metaES.se"] <- paste(abr,"metaES.se",sep="")
    names(ES.est)[names(ES.est)=="metaES.pval"] <- paste(abr,"metaES.pval",sep="")
    names(ES.est)[names(ES.est)=="metaES.ci.lb"] <- paste(abr,"metaES.ci.lb ",sep="")
    names(ES.est)[names(ES.est)=="metaES.ci.ub"] <- paste(abr,"metaES.ci.ub",sep="")
    
    cat("Recentered Effect Size Estimates", "\n")
    print(ES.est)
    return(list(rma.uncent=rma.uncent, ES.est=ES.est, 
                ES.mean=as.double(rma.noMed$b["intrcpt",1]), ES.SEM=as.double(rma.noMed$se[1])))
  }
}

#getmode
#Function to get the mode of a vector
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#IMPORT LAB DATA----

Lab <- read_excel("C:/Users/sbuja/Documents/Manuscripts for Publication/Quant Review Trans Addict Med/Meta Analysis/Lab Meta-Analysis Scored Data.xlsx", 
                 sheet = "Study Data", na="NA")

#How many medications?
length(table(Lab$Med))
#24 meds
#How many total effects
sum(table(Lab$Med))
#485

#CLEAN LAB DATA----

#Remove Self-Adminstration data
Lab.noSA <- subset(Lab, Admin != "Self-Admin")

#How many medications? 24
length(table(Lab.noSA$Med))
#How many total effects - 451 - 34 effects lost
sum(table(Lab.noSA$Med))

table(Lab.noSA$Admin)
#Challenge   Priming 
#      359        92

#exclude non main effects
Lab.noSA.main <- subset(Lab.noSA, PredTerm=="Med")

#How many medications? 24 - No medications lost
length(table(Lab.noSA.main$Med))
#How many total effects - 408 - 43 effects lost
sum(table(Lab.noSA.main$Med))

#How many studies? - 51
length(table(Lab.noSA.main$ID))

#How many samples? - 57
length(table(Lab.noSA.main$Sample))


#LAB STUDY STATS----
#Select sample level data for summary tables
Lab.Samples <- distinct(Lab.noSA.main[c("ID", "Sample", "Author", "Year", "Journal", "Med", "ActPlac", "MaxDose", "WithinMed", "Pharma",
                                        "N", "Admin", "Pop", "DpM", "IndSample", "MaxAlcDose")])
dim(Lab.Samples)
str(Lab.Samples)

#Medication
length(table(Lab.Samples$Med))
#24 total medications

table(Lab.Samples$Med)
# Acamprosate Aripiprazole     Baclofen  Dutasteride  Finasteride   Gabapentin    Ibudilast     Idazoxan Indomethacin   Isoflavone   Ivermectin 
#           2            2            2            1            1            2            1            1            1            1            1 
# Mecamylamine    Memantine    Nalmefene   Naltrexone   Olanzapine  Ondansetron   Paroxetine   Quetiapine   Rimonabant   Ritanserin   Topiramate 
#            3            2            2           23            2            1            1            1            1            1            1 
# Varenicline   Zonisamide 
#           3            1 

#Year
table(Lab.Samples$Year)
SpDesc(Lab.Samples$Year)
Lab.Years.Hist <- SpHist(Lab.Samples$Year, bins=10)
Lab.Years.Hist
#ggsave(Lab.Years.Hist, filename="Lab.Years.Hist.png", width = 6, height = 5, dpi=300)

#Making categories
Lab.Samples <- Lab.Samples %>% mutate(YearBins = cut(Year, breaks=c(-Inf, 1994, 1999, 2004, 2009, 2014, Inf), 
                                     labels=c("Pre 1995", "1995-1999", "2000-2004", "2005-2009", "2010-2014", "Post 2014")))
table(Lab.Samples$YearBins)
# Pre 1995 1995-1999 2000-2004 2005-2009 2010-2014 Post 2014 
#        1         8        14        16        14         3

#Journal
table(Lab.Samples$Journal)

#Within Subjects Design
table(Lab.Samples$WithinMed)
table(Lab.Samples$WithinMed)/sum(table(Lab.Samples$WithinMed))
#  0   1 
# 22  25
#39% 61% 

#Sample Size
SpDesc(Lab.Samples$N)
#   nbr.val        min        max     median       mean    SE.mean        var    std.dev 
# 57.000000  10.000000  90.000000  25.000000  32.929825   2.894106 477.423559  21.850024 
Lab.N.Hist <- SpHist(Lab.Samples$N, bins=10)
Lab.N.Hist
#ggsave(Lab.N.Hist, filename="Lab.N.Hist.png", width = 6, height = 5, dpi=300)

Lab.Samples <- Lab.Samples %>% mutate(NBins = cut(N, breaks=c(-Inf, 19, 39, 59, 79, Inf), 
                                                   labels=c("1-19", "20-39", "40-59", "60-79", "80+")))
table(Lab.Samples$NBins)
# 1-19 20-39 40-59 60-79   80+ 
#   20    19     8     7     3  

table(Lab.Samples$Admin)
# Challenge   Priming 
#        44        13

table(Lab.Samples$ActPlac)
table(Lab.Samples$ActPlac)/57
# 0  1 
#55  2 
#96% 4%

#Pharma funding
table(Lab.Samples$Pharma)
#  0  1 
# 48  9

#Population
Lab.Samples$Pop <- factor(Lab.Samples$Pop, levels=c("Light", "Heavy", "AUD"))
table(Lab.Samples$Pop)
#Light Heavy   AUD 
#   24    19    14 


#Mean Drinks per Month
SpDesc(Lab.Samples$DpM)
# nbr.val         min         max      median        mean     SE.mean         var     std.dev 
# 57.000000    5.600000  336.600000   51.400000   89.271369    9.363888 4997.896522   70.695803 
Lab.DpM.Hist <- SpHist(Lab.Samples$DpM, bins=10)
Lab.DpM.Hist
#ggsave(Lab.DpM.Hist, filename="Lab.DpM.Hist.png", width = 6, height = 5, dpi=300)
#qq plot to look for normality of DpM
SpQQPlot(Lab.Samples$DpM)
shapiro.test(Lab.Samples$DpM) #p < 0.0001
SpQQPlot(log(Lab.Samples$DpM))
shapiro.test(log(Lab.Samples$DpM)) #p = 0.01994
#Use log DpM
Lab.Samples$logDpM <- log(Lab.Samples$DpM)
Lab.noSA.main$logDpM <- log(Lab.noSA.main$DpM)


#Max Alcohol Dose
SpDesc(Lab.Samples$MaxAlcDose)
# nbr.val          min          max       median         mean      SE.mean          var      std.dev 
# 5.700000e+01 1.000000e-02 1.150000e-01 6.000000e-02 6.032523e-02 3.354349e-03 6.413445e-04 2.532478e-02 
Lab.MaxAlcDose.Hist <- SpHist(Lab.Samples$MaxAlcDose, bins=10)
Lab.MaxAlcDose.Hist
#ggsave(Lab.MaxAlcDose.Hist, filename="Lab.MaxAlcDose.Hist.png", width = 6, height = 5, dpi=300)
#qq plot to look for normality of DpM
SpQQPlot(Lab.Samples$MaxAlcDose)
shapiro.test(Lab.Samples$MaxAlcDose) #p .06541



#Effect size summaries

#total number of effects: 408
length(Lab.noSA.main$ES)

#Number of effects with no stat reported:
table(Lab.noSA.main$NoStat)
table(Lab.noSA.main$NoStat)/408
#   0   1 
# 240 168 
# 59% 41%

#Are certain meds more common with no Stat
table(Lab.noSA.main$Med, Lab.noSA.main$NoStat)
#Paroxetine, rimonabant, zonisamide, Indomethacin have no reported stats!!
Lab.noSA.main$NoStat.str <- ifelse(Lab.noSA.main$NoStat==0, "Stat Reported", "Stat Omitted")

#Outcome domains
table(Lab.noSA.main$OutDomain)
table(Lab.noSA.main$OutDomain)/408
# Craving     NegMood    Sedation Stimulation 
#      72          46         171         119 
#     18%         11%         42%         29%

OutDomain.plot <- ggplot(Lab.noSA.main, aes(OutDomain)) + geom_bar(aes(fill=NoStat.str), width = 0.8) + 
  ggtitle("Number of Effect Sizes Per Domain") + 
  scale_x_discrete("Outcome Domain") +
  SpTheme(legend.position = "right") + theme(legend.title = element_blank())
OutDomain.plot
#ggsave(OutDomain.plot, filename="OutDomain.plot.png", width = 6, height = 5, dpi = 300)

Lab.noSA.main$Med <- factor(Lab.noSA.main$Med)
Medication.NOutcomes.Plot <- ggplot(Lab.noSA.main, aes(Med)) + geom_bar(aes(fill=NoStat.str), width=0.8) +
  ggtitle("Number of Effect Sizes Per Medication") + 
  scale_x_discrete("Medication", limits = rev(levels(Lab.noSA.main$Med))) +
  coord_flip() + 
  SpTheme(legend.position = "right") + theme(legend.title = element_blank())
Medication.NOutcomes.Plot
#ggsave(Medication.NOutcomes.Plot, filename="Medication.NOutcomes.Plot.png", width = 5, height = 9, dpi = 300)


#META-ANALYSIS OF LABORATORY OUTCOMES----

#Center Covariates
#Admin - dummy code center at challenge
Lab.noSA.main$Admin.C <- ifelse(Lab.noSA.main$Admin=="Challenge", 0, 1)
table(Lab.noSA.main$Admin, Lab.noSA.main$Admin.C)
Lab.Samples$Admin.C <- ifelse(Lab.Samples$Admin=="Challenge", 0, 1)
table(Lab.Samples$Admin, Lab.Samples$Admin.C)

#MaxAlcDose - center at 0.06
Lab.noSA.main$MaxAlcDose.C <- Lab.noSA.main$MaxAlcDose - 0.06
SpDesc(Lab.noSA.main[c("MaxAlcDose", "MaxAlcDose.C")])
Lab.Samples$MaxAlcDose.C <- Lab.Samples$MaxAlcDose - 0.06
SpDesc(Lab.Samples[c("MaxAlcDose", "MaxAlcDose.C")])

#MaxDose - logbase2 center at modal dose of that medication (e.g. med with modal 50mg, 100 = 2, 25 = 0.5)
Lab.noSA.main$MaxDose.C <- NA
for(i in 1:dim(Lab.noSA.main)[1]){
  Lab.noSA.main$MaxDose.C[i] <- log(Lab.noSA.main$MaxDose[i] / getmode(subset(Lab.noSA.main,Med==Lab.noSA.main$Med[i])$MaxDose), 2)
}
SpHist(Lab.noSA.main$MaxDose.C)
Lab.Samples$MaxDose.C <- NA
for(i in 1:dim(Lab.Samples)[1]){
  Lab.Samples$MaxDose.C[i] <- log(Lab.Samples$MaxDose[i] / getmode(subset(Lab.Samples,Med==Lab.Samples$Med[i])$MaxDose), 2)
}
SpHist(Lab.Samples$MaxDose.C)


#CONSERVATIVE APROACH (No Stat = 0 effect size)----
#CRAVING OUTCOME - Conservative Approach (no stat = 0)----

#Subset Craving Outcomes
Lab.Craving <- subset(Lab.noSA.main, OutDomain=="Craving")
#reset levels of Med and Sample
Lab.Craving$Med <- factor(Lab.Craving$Med)
Lab.Craving$Sample <- factor(Lab.Craving$Sample)

#checks
dim(Lab.Craving) #72 effect sizes
table(Lab.Craving$Med)
length(table(Lab.Craving$Med)) #18 medications with craving outcomes
table(Lab.Craving$Sample) #which samples gave data
#Number of samples with craving data
length(levels(Lab.Craving$Sample)) # 43 samples with craving outcomes


#Aggregate craving effect sizes
Lab.Craving.ES <- agg(data=Lab.Craving, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(Lab.Craving.ES)[names(Lab.Craving.ES)=="id"] <- "Sample"
dim(Lab.Craving.ES)

#merge aggregated effect sizes with 
dim(Lab.Craving.ES)
dim(Lab.Samples)
Lab.Craving.ES <- inner_join(Lab.Craving.ES, Lab.Samples, by="Sample")
dim(Lab.Craving.ES) #43 effect sizes across samples

#reset levels of Med and Sample
Lab.Craving.ES$Med <- factor(Lab.Craving.ES$Med)
Lab.Craving.ES$Sample <- factor(Lab.Craving.ES$Sample)

#count number of Craving outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(Lab.Craving.ES)[1])
{
  Lab.Craving.ES$Outcomes[i] <- nrow(subset(Lab.Craving, Sample==Lab.Craving.ES$Sample[i]))
}


#Craving - RMA Analyses
rma.Craving<- rma.recenterMed.Lab(Lab.Craving.ES, abr="Cr.")

#Craving - Forest Plot
#Saving Size 8x7
forest.rma(rma.Craving$rma.uncent,
           slab = paste(Lab.Craving.ES$Author, Lab.Craving.ES$Year,sep=", "),
           ilab = cbind(as.character(Lab.Craving.ES$Med), Lab.Craving.ES$MaxDose, round(Lab.Craving.ES$DpM, 1), round(Lab.Craving.ES$MaxAlcDose, 3)),
           ilab.xpos = c(-3.3, -2.4, 1.2, 1.8),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.Craving.ES$Med),
           xlab="Hedge's G")
text(-4.9, 45, "Author(s) and Year", pos = 4, cex=.6)
text(-3.3, 45, "Medication", pos = 4, cex=.6)
text(-2.4, 45, "Dose", pos = 4, cex=.6)
text(1.2, 45, "DpM", pos = 4, cex=.6)
text(1.8, 45, "BrAC", pos = 4, cex=.6)
text(2.6, 45, "Hedge's G [95% CI]", pos = 4, cex=.6)
text(0, 47, "Alcohol Craving")

#funnel plot
Craving.Funnel <- gg.funnel(es=Lab.Craving.ES$es, es.var=Lab.Craving.ES$var, 
          mean.effect=rma.Craving$ES.mean, se.effect=rma.Craving$ES.SEM,
          title="Lab Outcomes - Alcohol Craving", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
          lab=factor(Lab.Craving.ES$Med), labsTitle="Medication")
Craving.Funnel
#ggsave(Craving.Funnel, filename="Craving.Funnel.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Craving$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = -0.9612, p = 0.3364

#Plot meta-analyzed effect sizes
rma.Craving$ES.est$Med <- factor(rma.Craving$ES.est$Med)
Craving.ES.Plot <- ggplot(rma.Craving$ES.est, aes(x=Cr.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = Cr.metaES - Cr.metaES.se, xmax = Cr.metaES + Cr.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Craving Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Craving$ES.est$Med))) +
  ggtitle("Craving Effect Sizes") +
  SpTheme()
Craving.ES.Plot
#ggsave(Craving.ES.Plot, filename="Craving.ES.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES <- rma.Craving$ES.est




#STIMULATION OUTCOME - Conservative Approach (no stat = 0)----

#Subset Stimulation Outcomes
Lab.Stimulation <- subset(Lab.noSA.main, OutDomain=="Stimulation")
#reset levels of Med and Sample
Lab.Stimulation$Med <- factor(Lab.Stimulation$Med)
Lab.Stimulation$Sample <- factor(Lab.Stimulation$Sample)

#checks
dim(Lab.Stimulation) #119 effect sizes
table(Lab.Stimulation$Med)
length(table(Lab.Stimulation$Med)) #24 medications with Stimulation outcomes
table(Lab.Stimulation$Sample) #which samples gave data
#Number of samples with Stimulation data
length(levels(Lab.Stimulation$Sample)) # 50 samples with Stimulation outcomes

#Aggregate Stimulation effect sizes
Lab.Stimulation.ES <- agg(data=Lab.Stimulation, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(Lab.Stimulation.ES)[names(Lab.Stimulation.ES)=="id"] <- "Sample"
dim(Lab.Stimulation.ES)

#merge aggregated effect sizes with 
dim(Lab.Stimulation.ES)
dim(Lab.Samples)
Lab.Stimulation.ES <- inner_join(Lab.Stimulation.ES, Lab.Samples, by="Sample")
dim(Lab.Stimulation.ES) #50 effect sizes across samples

#reset levels of Med and Sample
Lab.Stimulation.ES$Med <- factor(Lab.Stimulation.ES$Med)
Lab.Stimulation.ES$Sample <- factor(Lab.Stimulation.ES$Sample)

#count number of Stimulation outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(Lab.Stimulation.ES)[1])
{
  Lab.Stimulation.ES$Outcomes[i] <- nrow(subset(Lab.Stimulation, Sample==Lab.Stimulation.ES$Sample[i]))
}


#Stimulation - RMA Analyses
rma.Stimulation<- rma.recenterMed.Lab(Lab.Stimulation.ES, abr="St.")

#Stimulation - Forest Plot
#Saving Size 7x9
forest.rma(rma.Stimulation$rma.uncent,
           slab = paste(Lab.Stimulation.ES$Author, Lab.Stimulation.ES$Year,sep=", "),
           ilab = cbind(as.character(Lab.Stimulation.ES$Med), Lab.Stimulation.ES$MaxDose, round(Lab.Stimulation.ES$DpM, 1), round(Lab.Stimulation.ES$MaxAlcDose, 3)),
           ilab.xpos = c(-4.5, -3, 3.8, 4.7),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.Stimulation.ES$Med),
           xlab="Hedge's G")
text(-7.3, 53, "Author(s) and Year", pos = 4, cex=.6)
text(-4.5, 53, "Medication", pos = 4, cex=.6)
text(-3, 53, "Dose", pos = 4, cex=.6)
text(3.8, 53, "DpM", pos = 4, cex=.6)
text(4.7, 53, "BrAC", pos = 4, cex=.6)
text(6, 53, "Hedge's G [95% CI]", pos = 4, cex=.6)
text(0, 54.5, "Alcohol Stimulation")

#funnel plot
Stimulation.Funnel <- gg.funnel(es=Lab.Stimulation.ES$es, es.var=Lab.Stimulation.ES$var, 
                            mean.effect=rma.Stimulation$ES.mean, se.effect=rma.Stimulation$ES.SEM,
                            title="Lab Outcomes - Alcohol Stimulation", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                            lab=factor(Lab.Stimulation.ES$Med), labsTitle="Medication")
Stimulation.Funnel <- Stimulation.Funnel + annotate("rect", xmin = rma.Stimulation$ES.mean+1.96*0.6, xmax = 1.9, ymin = 0, ymax = 0.6, alpha = .1, fill="black")
Stimulation.Funnel
#ggsave(Stimulation.Funnel, filename="Stimulation.Funnel.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Stimulation$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = -0.8532, p = 0.3935

#Plot meta-analyzed effect sizes
rma.Stimulation$ES.est$Med <- factor(rma.Stimulation$ES.est$Med)
Stimulation.ES.Plot <- ggplot(rma.Stimulation$ES.est, aes(x=St.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = St.metaES - St.metaES.se, xmax = St.metaES + St.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Stimulation Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Stimulation$ES.est$Med))) +
  ggtitle("Stimulation Effect Sizes") +
  SpTheme()
Stimulation.ES.Plot
#ggsave(Stimulation.ES.Plot, filename="Stimulation.ES.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES <- full_join(Full.ES, rma.Stimulation$ES.est, by="Med")



#SEDATION OUTCOME - Conservative Approach (no stat = 0)----

#Subset Sedation Outcomes
Lab.Sedation <- subset(Lab.noSA.main, OutDomain=="Sedation")
#reset levels of Med and Sample
Lab.Sedation$Med <- factor(Lab.Sedation$Med)
Lab.Sedation$Sample <- factor(Lab.Sedation$Sample)

#checks
dim(Lab.Sedation) #171 effect sizes
table(Lab.Sedation$Med)
length(table(Lab.Sedation$Med)) #23 medications with Sedation outcomes
table(Lab.Sedation$Sample) #which samples gave data
#Number of samples with Sedation data
length(levels(Lab.Sedation$Sample)) # 53 samples with Sedation outcomes

#Aggregate Sedation effect sizes
Lab.Sedation.ES <- agg(data=Lab.Sedation, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(Lab.Sedation.ES)[names(Lab.Sedation.ES)=="id"] <- "Sample"
dim(Lab.Sedation.ES)

#merge aggregated effect sizes with 
dim(Lab.Sedation.ES)
dim(Lab.Samples)
Lab.Sedation.ES <- inner_join(Lab.Sedation.ES, Lab.Samples, by="Sample")
dim(Lab.Sedation.ES) #53 effect sizes across samples

#reset levels of Med and Sample
Lab.Sedation.ES$Med <- factor(Lab.Sedation.ES$Med)
Lab.Sedation.ES$Sample <- factor(Lab.Sedation.ES$Sample)

#count number of Sedation outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(Lab.Sedation.ES)[1])
{
  Lab.Sedation.ES$Outcomes[i] <- nrow(subset(Lab.Sedation, Sample==Lab.Sedation.ES$Sample[i]))
}


#Sedation - RMA Analyses
rma.Sedation<- rma.recenterMed.Lab(Lab.Sedation.ES, abr="Se.")

#Sedation - Forest Plot
#Saving Size 8x8
forest.rma(rma.Sedation$rma.uncent,
           slab = paste(Lab.Sedation.ES$Author, Lab.Sedation.ES$Year,sep=", "),
           ilab = cbind(as.character(Lab.Sedation.ES$Med), Lab.Sedation.ES$MaxDose, round(Lab.Sedation.ES$DpM, 1), round(Lab.Sedation.ES$MaxAlcDose, 3)),
           ilab.xpos = c(-2.7, -1.7, 2, 2.7),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.Sedation.ES$Med),
           xlab="Hedge's G")
text(-4.25, 56, "Author(s) and Year", pos = 4, cex=.6)
text(-2.7, 56, "Medication", pos = 4, cex=.6)
text(-1.7, 56, "Dose", pos = 4, cex=.6)
text(2, 56, "DpM", pos = 4, cex=.6)
text(2.7, 56, "BrAC", pos = 4, cex=.6)
text(3.3, 56, "Hedge's G [95% CI]", pos = 4, cex=.6)
text(0, 57.5, "Alcohol Sedation")

#funnel plot
Sedation.Funnel <- gg.funnel(es=Lab.Sedation.ES$es, es.var=Lab.Sedation.ES$var, 
                                mean.effect=rma.Sedation$ES.mean, se.effect=rma.Sedation$ES.SEM,
                                title="Lab Outcomes - Alcohol Sedation", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                                lab=factor(Lab.Sedation.ES$Med), labsTitle="Medication")
Sedation.Funnel
#ggsave(Sedation.Funnel, filename="Sedation.Funnel.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Sedation$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = 0.5219, p = 0.6018


#Plot meta-analyzed effect sizes
rma.Sedation$ES.est$Med <- factor(rma.Sedation$ES.est$Med)
Sedation.ES.Plot <- ggplot(rma.Sedation$ES.est, aes(x=Se.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = Se.metaES - Se.metaES.se, xmax = Se.metaES + Se.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Sedation Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Sedation$ES.est$Med))) +
  ggtitle("Sedation Effect Sizes") +
  SpTheme()
Sedation.ES.Plot
#ggsave(Sedation.ES.Plot, filename="Sedation.ES.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES <- full_join(Full.ES, rma.Sedation$ES.est, by="Med")



#NEGMOOD OUTCOME - Conservative Approach (no stat = 0)----

#Subset NegMood Outcomes
Lab.NegMood <- subset(Lab.noSA.main, OutDomain=="NegMood")
#reset levels of Med and Sample
Lab.NegMood$Med <- factor(Lab.NegMood$Med)
Lab.NegMood$Sample <- factor(Lab.NegMood$Sample)

#checks
dim(Lab.NegMood) #46 effect sizes
table(Lab.NegMood$Med)
length(table(Lab.NegMood$Med)) #12 medications with NegMood outcomes
table(Lab.NegMood$Sample) #which samples gave data
#Number of samples with NegMood data
length(levels(Lab.NegMood$Sample)) # 21 samples with NegMood outcomes


#Aggregate NegMood effect sizes
Lab.NegMood.ES <- agg(data=Lab.NegMood, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(Lab.NegMood.ES)[names(Lab.NegMood.ES)=="id"] <- "Sample"
dim(Lab.NegMood.ES)

#merge aggregated effect sizes with 
dim(Lab.NegMood.ES)
dim(Lab.Samples)
Lab.NegMood.ES <- inner_join(Lab.NegMood.ES, Lab.Samples, by="Sample")
dim(Lab.NegMood.ES) #21 effect sizes across samples

#reset levels of Med and Sample
Lab.NegMood.ES$Med <- factor(Lab.NegMood.ES$Med)
Lab.NegMood.ES$Sample <- factor(Lab.NegMood.ES$Sample)

#count number of NegMood outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(Lab.NegMood.ES)[1])
{
  Lab.NegMood.ES$Outcomes[i] <- nrow(subset(Lab.NegMood, Sample==Lab.NegMood.ES$Sample[i]))
}


#NegMood - RMA Analyses
rma.NegMood<- rma.recenterMed.Lab(Lab.NegMood.ES, abr="NM.")

#NegMood - Forest Plot
#Saving Size 8x7
forest.rma(rma.NegMood$rma.uncent,
           slab = paste(Lab.NegMood.ES$Author, Lab.NegMood.ES$Year,sep=", "),
           ilab = cbind(as.character(Lab.NegMood.ES$Med), Lab.NegMood.ES$MaxDose, round(Lab.NegMood.ES$DpM, 1), round(Lab.NegMood.ES$MaxAlcDose, 3)),
           ilab.xpos = c(-2.2, -1.2, 2, 2.7),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.NegMood.ES$Med),
           xlab="Hedge's G")
text(-4, 23, "Author(s) and Year", pos = 4, cex=1)
text(-2.2, 23, "Medication", pos = 4, cex=1)
text(-1.2, 23, "Dose", pos = 4, cex=1)
text(2, 23, "DpM", pos = 4, cex=1)
text(2.7, 23, "BrAC", pos = 4, cex=1)
text(3.3, 23, "Hedge's G [95% CI]", pos = 4, cex=1)
text(0, 24.5, "Alcohol Negative Mood", cex=1.1)

#funnel plot
NegMood.Funnel <- gg.funnel(es=Lab.NegMood.ES$es, es.var=Lab.NegMood.ES$var, 
                             mean.effect=rma.NegMood$ES.mean, se.effect=rma.NegMood$ES.SEM,
                             title="Lab Outcomes - Alcohol NegMood", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                             lab=factor(Lab.NegMood.ES$Med), labsTitle="Medication")
NegMood.Funnel <- NegMood.Funnel+ annotate("rect", xmin = rma.NegMood$ES.mean+1.96*0.4, xmax = 1.3, ymin = 0, ymax = 0.4, alpha = .1, fill="black")
NegMood.Funnel
#ggsave(NegMood.Funnel, filename="NegMood.Funnel.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.NegMood$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = 2.5325, p = 0.0113 ###ASYMMETRY TEST IS SIGNIFICANT###

#Plot meta-analyzed effect sizes
rma.NegMood$ES.est$Med <- factor(rma.NegMood$ES.est$Med)
NegMood.ES.Plot <- ggplot(rma.NegMood$ES.est, aes(x=NM.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = NM.metaES - NM.metaES.se, xmax = NM.metaES + NM.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("NegMood Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.NegMood$ES.est$Med))) +
  ggtitle("NegMood Effect Sizes") +
  SpTheme()
NegMood.ES.Plot
#ggsave(NegMood.ES.Plot, filename="NegMood.ES.Plot.png", width = 6, height = 5, dpi = 400)


#Save Medication Values
Full.ES <- full_join(Full.ES, rma.NegMood$ES.est, by="Med")

#Remaking Full.ES$Med a factor
Full.ES$Med <- factor(Full.ES$Med)




#IMPORT RCT DATA----
RCT <- read_excel("C:/Users/sbuja/Documents/Manuscripts for Publication/Quant Review Trans Addict Med/Meta Analysis/RCT Meta-Analysis Scored Data.xlsx", 
                  sheet = "Study Data", na="NA")

#How many medications?
length(table(RCT$Med))
#19 meds
#How many total effects
sum(table(RCT$Med))
#792


#CLEAN RCT DATA----

#Mean Impute DpM
mean(RCT$DpM, na.rm=T)
#266.4015

RCT$DpM <-ifelse(is.na(RCT$DpM), 266.4015, RCT$DpM)

str(RCT)

#delete rows from cochrane reviewed studies where the review did not include the outcome of interest
RCT.Clean <- subset(RCT, !is.na(NoStat))
dim(RCT.Clean) #663 effect sizes

#RCT STUDY STATS----
#Select sample level data for summary tables
RCT.Samples <- distinct(RCT[c("ID", "Cochrane", "Sample", "Author", "Year", "Journal", "Med", "ActPlac", "MaxDose", "Pharma",
                                        "N", "DpM", "Location", "AbsReq", "DepReq", "TrxDur", "FUDur")])

dim(RCT.Samples) #132 samples
str(RCT.Samples)

#Medication
length(table(RCT.Samples$Med))
#19 total medications

table(RCT.Samples$Med)
# Acamprosate  Aripiprazole      Baclofen Carbamazepine    Gabapentin Levetiracetam     Memantine     Nalmefene    Naltrexone    Olanzapine 
#          28             1             7             1             6             3             1             7            44             2 
# Ondansetron    Quetiapine    Rimonabant    Ritanserin    Sertraline    Topiramate     Valproate   Varenicline    Zonisamide 
#           3             5             1             3             1            10             3             4             2 

table(RCT.Samples$Cochrane)
# 0  1 
#56 76


#Year
table(RCT.Samples$Year)
SpDesc(RCT.Samples$Year)
# nbr.val          min          max       median         mean      SE.mean          var      std.dev 
# 132.0000000 1985.0000000 2016.0000000 2006.0000000 2005.5151515    0.5798973   44.3890817    6.6625132 
RCT.Years.Hist <- SpHist(RCT.Samples$Year, bins=10)
RCT.Years.Hist
#ggsave(RCT.Years.Hist, filename="RCT.Years.Hist.png", width = 6, height = 5, dpi=300)

#Making categories
RCT.Samples <- RCT.Samples %>% mutate(YearBins = cut(Year, breaks=c(-Inf, 1989, 1994, 1999, 2004, 2009, 2014, Inf), 
                                                     labels=c("Pre 1990", "1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014", "Post 2014")))
table(RCT.Samples$YearBins)
# Pre 1990 1990-1994 1995-1999 2000-2004 2005-2009 2010-2014 Post 2014 
#        1         5        19        33        33        30        11 

#ActPlac
table(RCT.Samples$ActPlac)
#  0   1 
#129   3 

#Pharma
table(RCT.Samples$Pharma)
# 0  1 
#58 59 

#N
table(RCT.Samples$N)
SpDesc(RCT.Samples$N)
#    nbr.val         min         max      median        mean     SE.mean         var     std.dev 
#  132.00000    10.00000  1383.00000   119.50000   183.97727    18.32964 44348.77047   210.59148
RCT.N.Hist <- SpHist(RCT.Samples$N, bins=20)
RCT.N.Hist + scale_x_continuous("N", breaks=seq(0,1400,200))
#ggsave(RCT.N.Hist, filename="RCT.N.Hist.png", width = 6, height = 5, dpi=300)
#Making categories
RCT.Samples <- RCT.Samples %>% mutate(NBins = cut(N, breaks=c(-Inf, 99, 199, 299, 399, 499, 599, Inf), 
                                                     labels=c("< 100", "100-199", "200-299", "300-399", "400-499", "500-599", ">600")))
table(RCT.Samples$NBins)
#  < 100 100-199 200-299 300-399 400-499 500-599    >600 
#     55      37      21       5       5       3       6 

#DpM
SpDesc(RCT.Samples$DpM)
#    nbr.val         min         max      median        mean     SE.mean         var     std.dev 
# 132.000000   68.600000  771.400000  266.401500  266.401511    8.599374 9761.298130   98.799282
RCT.DpM.Hist <- SpHist(RCT.Samples$DpM, bins=15)
RCT.DpM.Hist
#ggsave(RCT.DpM.Hist, filename="RCT.DpM.Hist.png", width = 6, height = 5, dpi=300)

#AbsReq
table(RCT.Samples$AbsReq)
# 0  1 
#81 36

#DepReq
table(RCT.Samples$DepReq)
#  0   1 
#  8 123 

#TrxDur
SpDesc(RCT.Samples$TrxDur)
#    nbr.val         min         max      median        mean     SE.mean         var     std.dev 
#132.0000000   4.0000000  52.0000000  12.0000000  15.5530303   0.8290299  90.7223572   9.5248285
RCT.TrxDur.Hist <- SpHist(RCT.Samples$TrxDur, bins=15)
RCT.TrxDur.Hist
#ggsave(RCT.TrxDur.Hist, filename="RCT.TrxDur.Hist.png", width = 6, height = 5, dpi=300)
RCT.Samples <- RCT.Samples %>% mutate(TrxDurBins = cut(TrxDur, breaks=c(-Inf, 6, 10, 14, 18, 22, 26, Inf), 
                                                  labels=c("<7", "7-10", "11-14", "15-18", "19-22", "23-26", ">26")))
table(RCT.Samples$TrxDurBins)
#   <7  7-10 11-14 15-18 19-22 23-26   >26 
#    6     4    90     6     1    18     7

#FUDur
SpDesc(RCT.Samples$FUDur)
#   nbr.val        min        max     median       mean    SE.mean        var    std.dev 
#123.000000   7.000000 104.000000  12.000000  22.691057   1.786787 392.690657  19.816424 
RCT.FUDur.Hist <- SpHist(RCT.Samples$FUDur, bins=15)
RCT.FUDur.Hist
#ggsave(RCT.FUDur.Hist, filename="RCT.FUDur.Hist.png", width = 6, height = 5, dpi=300)
RCT.Samples <- RCT.Samples %>% mutate(FUDurBins = cut(FUDur, breaks=c(-Inf, 12, 24, 36, 48, 60, 72, Inf)))
table(RCT.Samples$FUDurBins)
#(-Inf,12]   (12,24]   (24,36]   (36,48]   (48,60]   (60,72] (72, Inf] 
#       66        23         9         2        12         2         4


#Effect size summaries

#total number of effects: 663
length(RCT.Clean$ES)

#Number of effects with no stat reported:
table(RCT.Clean$NoStat)
table(RCT.Clean$NoStat)/664
#   0   1 
# 336 327 
# 51% 49% 

#How many meds have each outcome?
table(RCT.Clean$Med, RCT.Clean$OutName)
table(RCT.Clean$Med, RCT.Clean$OutDomain)

RCT.Clean$Cochrane.Str <- ifelse(RCT.Clean$Cochrane==1, "Cochrane", "Non-Cochrane")
RCT.Clean$Med <- factor(RCT.Clean$Med)
Medication.Cochrane.Plot <- ggplot(RCT.Clean, aes(Med)) + geom_bar(aes(fill=Cochrane.Str), width=0.8) +
  ggtitle("Number of RCT Effect Sizes Per Medication") + 
  scale_x_discrete("Medication", limits = rev(levels(RCT.Clean$Med))) +
  coord_flip() + 
  SpTheme(legend.position = "right") + theme(legend.title = element_blank())
Medication.Cochrane.Plot
#ggsave(Medication.Cochrane.Plot, filename="Medication.Cochrane.Plot.png", width = 5, height = 9, dpi = 300)


#META-ANALYSIS OF RCT OUTCOMES----

#Center Covariates - MaxDose.C, TrxDur.C, Pharma

#MaxDose - logbase2 center at modal dose of that medication (e.g. med with modal 50mg, 100 = 2, 25 = 0.5)
RCT.Clean$MaxDose.C <- NA
for(i in 1:dim(RCT.Clean)[1]){
  RCT.Clean$MaxDose.C[i] <- log(RCT.Clean$MaxDose[i] / getmode(subset(RCT.Clean,Med==RCT.Clean$Med[i])$MaxDose), 2)
}
SpHist(RCT.Clean$MaxDose.C)
RCT.Samples$MaxDose.C <- NA
for(i in 1:dim(RCT.Samples)[1]){
  RCT.Samples$MaxDose.C[i] <- log(RCT.Samples$MaxDose[i] / getmode(subset(RCT.Samples,Med==RCT.Samples$Med[i])$MaxDose), 2)
}
SpHist(RCT.Samples$MaxDose.C)

#TrxDur - treatment duration in weeks - centered at 12wks
RCT.Clean$TrxDur.C <- RCT.Clean$TrxDur - 12
RCT.Samples$TrxDur.C <- RCT.Samples$TrxDur - 12



#HEAVY DRINKING OUTCOME - Conservative Approach (No stat = 0)----

#subset Heavy Drinking outcomes
RCT.Heavy <- subset(RCT.Clean, OutDomain == "Heavy Drinking")
#reset levels of Med and Sample
RCT.Heavy$Med <- factor(RCT.Heavy$Med)
RCT.Heavy$Sample <- factor(RCT.Heavy$Sample)

#checks
dim(RCT.Heavy) #406 effect sizes
table(RCT.Heavy$Med)
length(table(RCT.Heavy$Med)) #19 medications with Heavy Drinking outcomes
table(RCT.Heavy$Sample) #which samples gave data
#Number of samples with NegMood data
length(levels(RCT.Heavy$Sample)) # 132 samples with Heavy Drinking outcomes

#Aggregate Heavy effect sizes
RCT.Heavy.ES <- agg(data=RCT.Heavy, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(RCT.Heavy.ES)[names(RCT.Heavy.ES)=="id"] <- "Sample"
dim(RCT.Heavy.ES) #132 effect sizes

#merge aggregated effect sizes with 
dim(RCT.Heavy.ES)
dim(Lab.Samples)
RCT.Heavy.ES <- inner_join(RCT.Heavy.ES, RCT.Samples, by="Sample")
dim(RCT.Heavy.ES) #132 effect sizes across samples

#reset levels of Med and Sample
RCT.Heavy.ES$Med <- factor(RCT.Heavy.ES$Med)
RCT.Heavy.ES$Sample <- factor(RCT.Heavy.ES$Sample)

#count number of Heavy outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(RCT.Heavy.ES)[1])
{
  RCT.Heavy.ES$Outcomes[i] <- nrow(subset(RCT.Heavy, Sample==RCT.Heavy.ES$Sample[i]))
}

#Heavy - RMA Analyses
rma.Heavy<- rma.recenterMed.RCT(RCT.Heavy.ES, abr="He.")

#Heavy - Forest Plot
#Saving Size 8x7
forest.rma(rma.Heavy$rma.uncent,
           slab = paste(RCT.Heavy.ES$Author, RCT.Heavy.ES$Year,sep=", "),
           ilab = cbind(as.character(RCT.Heavy.ES$Med), RCT.Heavy.ES$MaxDose, round(RCT.Heavy.ES$TrxDur, 1)),
           ilab.xpos = c(-4, -3, 2.5),
           ilab.pos=c(4,4,4),
           order=order(RCT.Heavy.ES$Med),
           xlab="Hedge's G", 
           cex=0.5)
text(-5.5, 135, "Author(s) and Year", pos = 4, cex=0.6)
text(-4, 135, "Medication", pos = 4, cex=0.6)
text(-3, 135, "Dose", pos = 4, cex=0.6)
text(2.5, 135, "Treatment Duration", pos = 4, cex=0.6)
text(4, 135, "Hedge's G [95% CI]", pos = 4, cex=0.6)
text(0, 137, "Heavy Drinking", cex=.8)

#funnel plot
Heavy.Funnel <- gg.funnel(es=RCT.Heavy.ES$es, es.var=RCT.Heavy.ES$var, 
                            mean.effect=rma.Heavy$ES.mean, se.effect=rma.Heavy$ES.SEM,
                            title="RCT Outcomes - Heavy Drinking", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                            lab=factor(RCT.Heavy.ES$Med), labsTitle="Medication")
Heavy.Funnel

#ggsave(Heavy.Funnel, filename="Heavy.Funnel.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Heavy$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = -1.9100, p = 0.0561

#Plot meta-analyzed effect sizes
rma.Heavy$ES.est$Med <- factor(rma.Heavy$ES.est$Med)
Heavy.ES.Plot <- ggplot(rma.Heavy$ES.est, aes(x=He.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = He.metaES - He.metaES.se, xmax = He.metaES + He.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Heavy Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Heavy$ES.est$Med))) +
  ggtitle("Heavy Effect Sizes") +
  SpTheme()
Heavy.ES.Plot
#ggsave(Heavy.ES.Plot, filename="Heavy.ES.Plot.png", width = 6, height = 5, dpi = 400)


#Save Medication Values
Full.ES <- full_join(Full.ES, rma.Heavy$ES.est, by="Med")


#ABSTINENCE OUTCOME - Conservative Approach (No stat = 0)----

#subset Abstinence outcomes
RCT.Abstinence <- subset(RCT.Clean, OutDomain == "Abstinence")
#reset levels of Med and Sample
RCT.Abstinence$Med <- factor(RCT.Abstinence$Med)
RCT.Abstinence$Sample <- factor(RCT.Abstinence$Sample)

#checks
dim(RCT.Abstinence) #257 effect sizes
table(RCT.Abstinence$Med)
length(table(RCT.Abstinence$Med)) #19 medications with Abstinence outcomes
table(RCT.Abstinence$Sample) #which samples gave data
#Number of samples with NegMood data
length(levels(RCT.Abstinence$Sample)) # 129 samples with Abstinence outcomes

#Aggregate Abstinence effect sizes
RCT.Abstinence.ES <- agg(data=RCT.Abstinence, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(RCT.Abstinence.ES)[names(RCT.Abstinence.ES)=="id"] <- "Sample"
dim(RCT.Abstinence.ES)

#merge aggregated effect sizes with 
dim(RCT.Abstinence.ES)
dim(Lab.Samples)
RCT.Abstinence.ES <- inner_join(RCT.Abstinence.ES, RCT.Samples, by="Sample")
dim(RCT.Abstinence.ES) #129 effect sizes across samples

#reset levels of Med and Sample
RCT.Abstinence.ES$Med <- factor(RCT.Abstinence.ES$Med)
RCT.Abstinence.ES$Sample <- factor(RCT.Abstinence.ES$Sample)

#count number of Abstinence outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(RCT.Abstinence.ES)[1])
{
  RCT.Abstinence.ES$Outcomes[i] <- nrow(subset(RCT.Abstinence, Sample==RCT.Abstinence.ES$Sample[i]))
}

#Abstinence - RMA Analyses
rma.Abstinence<- rma.recenterMed.RCT(RCT.Abstinence.ES, abr="Ab.")

#Abstinence - Forest Plot
#Saving Size 8x7
forest.rma(rma.Abstinence$rma.uncent,
           slab = paste(RCT.Abstinence.ES$Author, RCT.Abstinence.ES$Year,sep=", "),
           ilab = cbind(as.character(RCT.Abstinence.ES$Med), RCT.Abstinence.ES$MaxDose, round(RCT.Abstinence.ES$TrxDur, 1)),
           ilab.xpos = c(-4, -2.5, 5),
           ilab.pos=c(4,4,4),
           order=order(RCT.Abstinence.ES$Med),
           xlab="Hedge's G", 
           cex=0.5)
text(-7.5, 131, "Author(s) and Year", pos = 4, cex=0.6)
text(-4, 131, "Medication", pos = 4, cex=0.6)
text(-2.5, 131, "Dose", pos = 4, cex=0.6)
text(5, 131, "Treatment Duration", pos = 4, cex=0.6)
text(8, 131, "Hedge's G [95% CI]", pos = 4, cex=0.6)
text(0, 133, "Abstinence", cex=.8)

#funnel plot
Abstinence.Funnel <- gg.funnel(es=RCT.Abstinence.ES$es, es.var=RCT.Abstinence.ES$var, 
                          mean.effect=rma.Abstinence$ES.mean, se.effect=rma.Abstinence$ES.SEM,
                          title="RCT Outcomes - Abstinence", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                          lab=factor(RCT.Abstinence.ES$Med), labsTitle="Medication")
Abstinence.Funnel + annotate("rect", xmin = rma.Abstinence$ES.mean+1.96*0.6, xmax = 3.2, ymin = 0, ymax = 0.6, alpha = .1, fill="black")

#ggsave(Abstinence.Funnel, filename="Abstinence.Funnel.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Abstinence$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = 2.0116, p = 0.0443 #Significant assymmetry on Abstinence outcomes. 

#Plot meta-analyzed effect sizes
rma.Abstinence$ES.est$Med <- factor(rma.Abstinence$ES.est$Med)
Abstinence.ES.Plot <- ggplot(rma.Abstinence$ES.est, aes(x=Ab.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = Ab.metaES - Ab.metaES.se, xmax = Ab.metaES + Ab.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Abstinence Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Abstinence$ES.est$Med))) +
  ggtitle("Abstinence Effect Sizes") +
  SpTheme()
Abstinence.ES.Plot
#ggsave(Abstinence.ES.Plot, filename="Abstinence.ES.Plot.png", width = 6, height = 5, dpi = 400)


#Save Medication Values
Full.ES <- full_join(Full.ES, rma.Abstinence$ES.est, by="Med")


#WILLIAMSON-YORK TRANSLATIONAL ANALYSIS - Conservative Approach (No Stat = 0)----

#Lab Craving - Heavy Drinking----
Cr.He.ES <- na.exclude(Full.ES[,c("Med", "Cr.metaES", "Cr.metaES.se", "He.metaES", "He.metaES.se")])
dim(Cr.He.ES)
#13 medications for this analysis

Cr.He.WYbwls <- WYbwls(x=Cr.He.ES$Cr.metaES, xsd=Cr.He.ES$Cr.metaES.se,
                       y=Cr.He.ES$He.metaES, ysd=Cr.He.ES$He.metaES.se,
                       print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#             Est  	   SE 
# Int   	 -0.798 	 0.918 
# Slope 	 -3.6356 	 4.689 
# 
# r:   	 -0.1617 
# r^2: 	 0.0262 
# p:   	 0.454529412693122

#Modify plot with specific labels
Cr.He.WYbwls.plot <- Cr.He.WYbwls$plot + 
  ggtitle("Laboratory Craving and RCT Heavy Drinking\n- Conservative Approach -") +
  scale_x_continuous("Laboratory Effects on Alcohol Craving (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") +
  theme(legend.position=c(0.2,0.2))
Cr.He.WYbwls.plot
#ggsave(Cr.He.WYbwls.plot, filename="Cr.He.WYbwls.plot.png", width = 6, height = 5, dpi = 400)


#Lab Stimulation - Heavy Drinking----
St.He.ES <- na.exclude(Full.ES[,c("Med", "St.metaES", "St.metaES.se", "He.metaES", "He.metaES.se")])
dim(St.He.ES)
#15 medications for this analysis

St.He.WYbwls <- WYbwls(x=St.He.ES$St.metaES, xsd=St.He.ES$St.metaES.se,
                       y=St.He.ES$He.metaES, ysd=St.He.ES$He.metaES.se,
                       print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 -0.303 	 0.255 
# Slope 	 2.3831 	 2.51 
# 
# r:   	 0.128 
# r^2: 	 0.0164 
# p:   	 0.359618589952904 

#Modify plot with specific labels
St.He.WYbwls.plot <- St.He.WYbwls$plot + 
  ggtitle("Laboratory Stimulation and RCT Heavy Drinking\n- Conservative Approach -") +
  scale_x_continuous("Laboratory Effects on Alcohol Stimulation (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") +
  theme(legend.position=c(.8,.6))
St.He.WYbwls.plot
#ggsave(St.He.WYbwls.plot, filename="St.He.WYbwls.plot.png", width = 6, height = 5, dpi = 400)


#Lab Sedation - Heavy Drinking----
Se.He.ES <- na.exclude(Full.ES[,c("Med", "Se.metaES", "Se.metaES.se", "He.metaES", "He.metaES.se")])
dim(Se.He.ES)
#15 medications for this analysis

Se.He.WYbwls <- WYbwls(x=Se.He.ES$Se.metaES, xsd=Se.He.ES$Se.metaES.se,
                       y=Se.He.ES$He.metaES, ysd=Se.He.ES$He.metaES.se,
                       print=T, plot=T, tol=1e-08)

#Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 -0.081 	 0.234 
# Slope 	 -4.0355 	 5.545 
# 
# r:   	 -0.3215 
# r^2: 	 0.1033 
# p:   	 0.479688870163164 

#Modify plot with specific labels
Se.He.WYbwls.plot <- Se.He.WYbwls$plot + 
  ggtitle("Laboratory Sedation and RCT Heavy Drinking\n- Conservative Approach -") +
  scale_x_continuous("Laboratory Effects on Alcohol Sedation (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") +
  theme(legend.position=c(0.2,0.2))
Se.He.WYbwls.plot
#ggsave(Se.He.WYbwls.plot, filename="Se.He.WYbwls.plot.png", width = 6, height = 5, dpi = 400)


#Lab NegMood - Heavy Drinking----
NM.He.ES <- na.exclude(Full.ES[,c("Med", "NM.metaES", "NM.metaES.se", "He.metaES", "He.metaES.se")])
dim(NM.He.ES)
#8 medications for this analysis

NM.He.WYbwls <- WYbwls(x=NM.He.ES$NM.metaES, xsd=NM.He.ES$NM.metaES.se,
                       y=NM.He.ES$He.metaES, ysd=NM.He.ES$He.metaES.se,
                       print=T, plot=T, tol=1e-08)

#  Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 -45.348 	 10694.719 
# Slope 	 362.1938 	 85627.106 
# 
# r:   	 -0.0779 
# r^2: 	 0.0061 
# p:   	 0.996762171256329 

#Modify plot with specific labels
NM.He.WYbwls.plot <- NM.He.WYbwls$plot + 
  ggtitle("Laboratory Negative Mood and RCT Heavy Drinking\n- Conservative Approach -") +
  scale_x_continuous("Laboratory Effects on Negative Mood (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)")
NM.He.WYbwls.plot
#ggsave(NM.He.WYbwls.plot, filename="NM.He.WYbwls.plot.png", width = 6, height = 5, dpi = 400)


#Lab Craving - Abstinence----
Cr.Ab.ES <- na.exclude(Full.ES[,c("Med", "Cr.metaES", "Cr.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(Cr.Ab.ES)
#13 medications for this analysis

Cr.Ab.WYbwls <- WYbwls(x=Cr.Ab.ES$Cr.metaES, xsd=Cr.Ab.ES$Cr.metaES.se,
                       y=Cr.Ab.ES$Ab.metaES, ysd=Cr.Ab.ES$Ab.metaES.se,
                       print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
#  Coefficients: 
#           Est  	   SE 
# Int   	 0.049 	 0.052 
# Slope 	 -0.2931 	 0.208 
# 
# r:   	 0.3239 
# r^2: 	 0.1049 
# p:   	 0.186989042274241 

#Modify plot with specific labels
Cr.Ab.WYbwls.plot <- Cr.Ab.WYbwls$plot + 
  ggtitle("Laboratory Craving and RCT Abstinence\n- Conservative Approach -") +
  scale_x_continuous("Laboratory Effects on Alcohol Craving (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)") 
Cr.Ab.WYbwls.plot
#ggsave(Cr.Ab.WYbwls.plot, filename="Cr.Ab.WYbwls.plot.png", width = 6, height = 5, dpi = 400)


#Lab Stimulation - Abstinence----
St.Ab.ES <- na.exclude(Full.ES[,c("Med", "St.metaES", "St.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(St.Ab.ES)
#15 medications for this analysis

St.Ab.WYbwls <- WYbwls(x=St.Ab.ES$St.metaES, xsd=St.Ab.ES$St.metaES.se,
                       y=St.Ab.ES$Ab.metaES, ysd=St.Ab.ES$Ab.metaES.se,
                       print=T, plot=T, tol=1e-08)

#  Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 0.093 	 0.039 
# Slope 	 -0.2437 	 0.149 
# 
# r:   	 -0.1565 
# r^2: 	 0.0245 
# p:   	 0.12604055854209 

#Modify plot with specific labels
St.Ab.WYbwls.plot <- St.Ab.WYbwls$plot + 
  ggtitle("Laboratory Stimulation and RCT Abstinence\n- Conservative Approach -") +
  scale_x_continuous("Laboratory Effects on Alcohol Stimulation (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
St.Ab.WYbwls.plot
#ggsave(St.Ab.WYbwls.plot, filename="St.Ab.WYbwls.plot.png", width = 6, height = 5, dpi = 400)


#Lab Sedation - Abstinence----
Se.Ab.ES <- na.exclude(Full.ES[,c("Med", "Se.metaES", "Se.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(Se.Ab.ES)
#15 medications for this analysis

Se.Ab.WYbwls <- WYbwls(x=Se.Ab.ES$Se.metaES, xsd=Se.Ab.ES$Se.metaES.se,
                       y=Se.Ab.ES$Ab.metaES, ysd=Se.Ab.ES$Ab.metaES.se,
                       print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 0.075 	 0.046 
# Slope 	 0.4586 	 0.267 
# 
# r:   	 0.4463 
# r^2: 	 0.1992 
# p:   	 0.109057436537707 

#Modify plot with specific labels
Se.Ab.WYbwls.plot <- Se.Ab.WYbwls$plot + 
  ggtitle("Laboratory Sedation and RCT Abstinence\n- Conservative Approach -") +
  scale_x_continuous("Laboratory Effects on Alcohol Sedation (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
Se.Ab.WYbwls.plot
#ggsave(Se.Ab.WYbwls.plot, filename="Se.Ab.WYbwls.plot.png", width = 6, height = 5, dpi = 400)


#Lab NegMood - Abstinence----
NM.Ab.ES <- na.exclude(Full.ES[,c("Med", "NM.metaES", "NM.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(NM.Ab.ES)
#8 medications for this analysis

NM.Ab.WYbwls <- WYbwls(x=NM.Ab.ES$NM.metaES, xsd=NM.Ab.ES$NM.metaES.se,
                       y=NM.Ab.ES$Ab.metaES, ysd=NM.Ab.ES$Ab.metaES.se,
                       print=T, plot=T, tol=1e-08)

#  Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
#  Coefficients: 
# Est  	   SE 
# Int   	 -0.108 	 0.245 
# Slope 	 1.3301 	 1.483 
# 
# r:   	 0.5049 
# r^2: 	 0.2549 
# p:   	 0.404255122726759 

#Modify plot with specific labels
NM.Ab.WYbwls.plot <- NM.Ab.WYbwls$plot + 
  ggtitle("Laboratory Negative Mood and RCT Abstinence\n- Conservative Approach -") +
  scale_x_continuous("Laboratory Effects on Negative Mood (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
NM.Ab.WYbwls.plot
#ggsave(NM.Ab.WYbwls.plot, filename="NM.Ab.WYbwls.plot.png", width = 6, height = 5, dpi = 400)



#MODERATE APROACH (No Stat = effect size associated with p-value of 0.5)----

#CRAVING OUTCOME - Moderate Approach----

#impute no stat effect sizes to p=0.5
Lab.Craving$ES.Imp <- NA
Lab.Craving$ES.Impvar <- NA

for(i in 1:dim(Lab.Craving)[1]){
  if(Lab.Craving$NoStat[i]==1){
    Lab.Craving$ES.Imp[i] <- -as.double(pes(p=0.5, n.1=Lab.Craving$N[i]/2, n.2=Lab.Craving$N[i]/2, dig=4, verbose=F)[c("g")])
    Lab.Craving$ES.Impvar[i] <- as.double(pes(p=0.5, n.1=Lab.Craving$N[i]/2, n.2=Lab.Craving$N[i]/2, dig=4, verbose=F)[c("var.g")])
  }
  else{
    Lab.Craving$ES.Imp[i] <- Lab.Craving$ES[i]
    Lab.Craving$ES.Impvar[i] <- Lab.Craving$ESvar[i]
  }
}

#checks
length(Lab.Craving$ES.Imp) #72 effect sizes
length(table(Lab.Craving$Med)) #18 medications with craving outcomes
table(Lab.Craving$Sample) #which samples gave data
#Number of samples with craving data
length(levels(Lab.Craving$Sample)) # 43 samples with craving outcomes

#Aggregate craving effect sizes
Lab.Craving.ES.Imp <- agg(data=Lab.Craving, id=Sample, es=ES.Imp, var=ES.Impvar,  method = "BHHR", cor=.6)
names(Lab.Craving.ES.Imp)[names(Lab.Craving.ES.Imp)=="id"] <- "Sample"
dim(Lab.Craving.ES.Imp) #43 effect sizes

#merge aggregated effect sizes with 
dim(Lab.Craving.ES.Imp)
dim(Lab.Samples)
Lab.Craving.ES.Imp <- inner_join(Lab.Craving.ES.Imp, Lab.Samples, by="Sample")
dim(Lab.Craving.ES.Imp) #43 effect sizes across samples

#reset levels of Med and Sample
Lab.Craving.ES.Imp$Med <- factor(Lab.Craving.ES.Imp$Med)
Lab.Craving.ES.Imp$Sample <- factor(Lab.Craving.ES.Imp$Sample)

#count number of Craving outcomes that were aggregated for each aggregated ES.Imp
#count outcomes
for (i in 1:dim(Lab.Craving.ES.Imp)[1])
{
  Lab.Craving.ES.Imp$Outcomes[i] <- nrow(subset(Lab.Craving, Sample==Lab.Craving.ES.Imp$Sample[i]))
}


#Craving - RMA Analyses
rma.Craving.Imp<- rma.recenterMed.Lab(Lab.Craving.ES.Imp, abr="Cr.")

#Craving - Forest Plot
#Saving Size 8x7
forest.rma(rma.Craving.Imp$rma.uncent,
           slab = paste(Lab.Craving.ES.Imp$Author, Lab.Craving.ES.Imp$Year,sep=", "),
           ilab = cbind(as.character(Lab.Craving.ES.Imp$Med), Lab.Craving.ES.Imp$MaxDose, round(Lab.Craving.ES.Imp$DpM, 1), round(Lab.Craving.ES.Imp$MaxAlcDose, 3)),
           ilab.xpos = c(-3.3, -2.4, 1.2, 1.8),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.Craving.ES.Imp$Med),
           xlab="Hedge's G")
text(-4.9, 45, "Author(s) and Year", pos = 4, cex=.6)
text(-3.3, 45, "Medication", pos = 4, cex=.6)
text(-2.4, 45, "Dose", pos = 4, cex=.6)
text(1.2, 45, "DpM", pos = 4, cex=.6)
text(1.8, 45, "BrAC", pos = 4, cex=.6)
text(2.6, 45, "Hedge's G [95% CI]", pos = 4, cex=.6)
text(0, 47, "Alcohol Craving - Moderate Imputation")

#funnel plot
Craving.Funnel.Imp <- gg.funnel(es=Lab.Craving.ES.Imp$es, es.var=Lab.Craving.ES.Imp$var, 
                            mean.effect=rma.Craving.Imp$ES.mean, se.effect=rma.Craving.Imp$ES.SEM,
                            title="Lab Outcomes - Alcohol Craving\nModerate Imputation", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                            lab=factor(Lab.Craving.ES.Imp$Med), labsTitle="Medication")
Craving.Funnel.Imp
#ggsave(Craving.Funnel.Imp, filename="Craving.Funnel.Imp.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Craving.Imp$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = -0.3364, p = 0.7365

#Plot meta-analyzed effect sizes
rma.Craving.Imp$ES.est$Med <- factor(rma.Craving.Imp$ES.est$Med)
Craving.ES.Imp.Plot <- ggplot(rma.Craving.Imp$ES.est, aes(x=Cr.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = Cr.metaES - Cr.metaES.se, xmax = Cr.metaES + Cr.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Craving Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Craving.Imp$ES.est$Med))) +
  ggtitle("Craving Effect Sizes\nModerate Imputation") +
  SpTheme()
Craving.ES.Imp.Plot
#ggsave(Craving.ES.Imp.Plot, filename="Craving.ES.Imp.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES.Imp <- rma.Craving.Imp$ES.est


#STIMULATION OUTCOME - Moderate Approach----

#impute no stat effect sizes to p=0.5
Lab.Stimulation$ES.Imp <- NA
Lab.Stimulation$ES.Impvar <- NA

for(i in 1:dim(Lab.Stimulation)[1]){
  if(Lab.Stimulation$NoStat[i]==1){
    Lab.Stimulation$ES.Imp[i] <- -as.double(pes(p=0.5, n.1=Lab.Stimulation$N[i]/2, n.2=Lab.Stimulation$N[i]/2, dig=4, verbose=F)[c("g")])
    Lab.Stimulation$ES.Impvar[i] <- as.double(pes(p=0.5, n.1=Lab.Stimulation$N[i]/2, n.2=Lab.Stimulation$N[i]/2, dig=4, verbose=F)[c("var.g")])
  }
  else{
    Lab.Stimulation$ES.Imp[i] <- Lab.Stimulation$ES[i]
    Lab.Stimulation$ES.Impvar[i] <- Lab.Stimulation$ESvar[i]
  }
}

#checks
length(Lab.Stimulation$ES.Imp) #119 effect sizes
length(table(Lab.Stimulation$Med)) #24 medications with craving outcomes
table(Lab.Stimulation$Sample) #which samples gave data
#Number of samples with craving data
length(levels(Lab.Stimulation$Sample)) # 50 samples with craving outcomes

#Aggregate craving effect sizes
Lab.Stimulation.ES.Imp <- agg(data=Lab.Stimulation, id=Sample, es=ES.Imp, var=ES.Impvar,  method = "BHHR", cor=.6)
names(Lab.Stimulation.ES.Imp)[names(Lab.Stimulation.ES.Imp)=="id"] <- "Sample"
dim(Lab.Stimulation.ES.Imp) #50 effect sizes

#merge aggregated effect sizes with 
dim(Lab.Stimulation.ES.Imp)
dim(Lab.Samples)
Lab.Stimulation.ES.Imp <- inner_join(Lab.Stimulation.ES.Imp, Lab.Samples, by="Sample")
dim(Lab.Stimulation.ES.Imp) #50 effect sizes across samples

#reset levels of Med and Sample
Lab.Stimulation.ES.Imp$Med <- factor(Lab.Stimulation.ES.Imp$Med)
Lab.Stimulation.ES.Imp$Sample <- factor(Lab.Stimulation.ES.Imp$Sample)

#count number of Stimulation outcomes that were aggregated for each aggregated ES.Imp
#count outcomes
for (i in 1:dim(Lab.Stimulation.ES.Imp)[1])
{
  Lab.Stimulation.ES.Imp$Outcomes[i] <- nrow(subset(Lab.Stimulation, Sample==Lab.Stimulation.ES.Imp$Sample[i]))
}


#Stimulation - RMA Analyses
rma.Stimulation.Imp<- rma.recenterMed.Lab(Lab.Stimulation.ES.Imp, abr="St.")

#Stimulation - Forest Plot
#Saving Size 7x9
forest.rma(rma.Stimulation.Imp$rma.uncent,
           slab = paste(Lab.Stimulation.ES.Imp$Author, Lab.Stimulation.ES.Imp$Year,sep=", "),
           ilab = cbind(as.character(Lab.Stimulation.ES.Imp$Med), Lab.Stimulation.ES.Imp$MaxDose, round(Lab.Stimulation.ES.Imp$DpM, 1), round(Lab.Stimulation.ES.Imp$MaxAlcDose, 3)),
           ilab.xpos = c(-4.5, -3, 3.8, 4.7),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.Stimulation.ES.Imp$Med),
           xlab="Hedge's G")
text(-7.3, 53, "Author(s) and Year", pos = 4, cex=.6)
text(-4.5, 53, "Medication", pos = 4, cex=.6)
text(-3, 53, "Dose", pos = 4, cex=.6)
text(3.8, 53, "DpM", pos = 4, cex=.6)
text(4.7, 53, "BrAC", pos = 4, cex=.6)
text(6, 53, "Hedge's G [95% CI]", pos = 4, cex=.6)
text(0, 54.5, "Alcohol Stimulation - Moderate Imputation")

#funnel plot
Stimulation.Funnel.Imp <- gg.funnel(es=Lab.Stimulation.ES.Imp$es, es.var=Lab.Stimulation.ES.Imp$var, 
                            mean.effect=rma.Stimulation.Imp$ES.mean, se.effect=rma.Stimulation.Imp$ES.SEM,
                            title="Lab Outcomes - Alcohol Stimulation\nModerate Imputation", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                            lab=factor(Lab.Stimulation.ES.Imp$Med), labsTitle="Medication")
Stimulation.Funnel.Imp<- Stimulation.Funnel.Imp  + annotate("rect", xmin = rma.Stimulation.Imp$ES.mean+1.96*0.6, xmax = 1.9, ymin = 0, ymax = 0.6, alpha = .1, fill="black")
Stimulation.Funnel.Imp
#ggsave(Stimulation.Funnel.Imp, filename="Stimulation.Funnel.Imp.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Stimulation.Imp$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = -0.8947, p = 0.3709

#Plot meta-analyzed effect sizes
rma.Stimulation.Imp$ES.est$Med <- factor(rma.Stimulation.Imp$ES.est$Med)
Stimulation.ES.Imp.Plot <- ggplot(rma.Stimulation.Imp$ES.est, aes(x=St.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = St.metaES - St.metaES.se, xmax = St.metaES + St.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Stimulation Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Stimulation.Imp$ES.est$Med))) +
  ggtitle("Stimulation Effect Sizes\nModerate Imputation") +
  SpTheme()
Stimulation.ES.Imp.Plot
#ggsave(Stimulation.ES.Imp.Plot, filename="Stimulation.ES.Imp.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES.Imp <- full_join(Full.ES.Imp, rma.Stimulation.Imp$ES.est, by="Med")


#SEDATION OUTCOME - Moderate Approach----

#impute no stat effect sizes to p=0.5
Lab.Sedation$ES.Imp <- NA
Lab.Sedation$ES.Impvar <- NA

for(i in 1:dim(Lab.Sedation)[1]){
  if(Lab.Sedation$NoStat[i]==1){
    Lab.Sedation$ES.Imp[i] <- as.double(pes(p=0.5, n.1=Lab.Sedation$N[i]/2, n.2=Lab.Sedation$N[i]/2, dig=4, verbose=F)[c("g")])
    Lab.Sedation$ES.Impvar[i] <- as.double(pes(p=0.5, n.1=Lab.Sedation$N[i]/2, n.2=Lab.Sedation$N[i]/2, dig=4, verbose=F)[c("var.g")])
  }
  else{
    Lab.Sedation$ES.Imp[i] <- Lab.Sedation$ES[i]
    Lab.Sedation$ES.Impvar[i] <- Lab.Sedation$ESvar[i]
  }
}

#checks
length(Lab.Sedation$ES.Imp) #171 effect sizes
length(table(Lab.Sedation$Med)) #23 medications with craving outcomes
table(Lab.Sedation$Sample) #which samples gave data
#Number of samples with craving data
length(levels(Lab.Sedation$Sample)) # 53 samples with craving outcomes

#Aggregate craving effect sizes
Lab.Sedation.ES.Imp <- agg(data=Lab.Sedation, id=Sample, es=ES.Imp, var=ES.Impvar,  method = "BHHR", cor=.6)
names(Lab.Sedation.ES.Imp)[names(Lab.Sedation.ES.Imp)=="id"] <- "Sample"
dim(Lab.Sedation.ES.Imp) #53 effect sizes

#merge aggregated effect sizes with 
dim(Lab.Sedation.ES.Imp)
dim(Lab.Samples)
Lab.Sedation.ES.Imp <- inner_join(Lab.Sedation.ES.Imp, Lab.Samples, by="Sample")
dim(Lab.Sedation.ES.Imp) #53 effect sizes across samples

#reset levels of Med and Sample
Lab.Sedation.ES.Imp$Med <- factor(Lab.Sedation.ES.Imp$Med)
Lab.Sedation.ES.Imp$Sample <- factor(Lab.Sedation.ES.Imp$Sample)

#count number of Sedation outcomes that were aggregated for each aggregated ES.Imp
#count outcomes
for (i in 1:dim(Lab.Sedation.ES.Imp)[1])
{
  Lab.Sedation.ES.Imp$Outcomes[i] <- nrow(subset(Lab.Sedation, Sample==Lab.Sedation.ES.Imp$Sample[i]))
}


#Sedation - RMA Analyses
rma.Sedation.Imp<- rma.recenterMed.Lab(Lab.Sedation.ES.Imp, abr="Se.")

#Sedation - Forest Plot
#Saving Size 8x8
forest.rma(rma.Sedation.Imp$rma.uncent,
           slab = paste(Lab.Sedation.ES.Imp$Author, Lab.Sedation.ES.Imp$Year,sep=", "),
           ilab = cbind(as.character(Lab.Sedation.ES.Imp$Med), Lab.Sedation.ES.Imp$MaxDose, round(Lab.Sedation.ES.Imp$DpM, 1), round(Lab.Sedation.ES.Imp$MaxAlcDose, 3)),
           ilab.xpos = c(-2.7, -1.7, 2, 2.7),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.Sedation.ES.Imp$Med),
           xlab="Hedge's G")
text(-4.25, 56, "Author(s) and Year", pos = 4, cex=.6)
text(-2.7, 56, "Medication", pos = 4, cex=.6)
text(-1.7, 56, "Dose", pos = 4, cex=.6)
text(2, 56, "DpM", pos = 4, cex=.6)
text(2.7, 56, "BrAC", pos = 4, cex=.6)
text(3.3, 56, "Hedge's G [95% CI]", pos = 4, cex=.6)
text(0, 57.5, "Alcohol Sedation - Moderate Imputation")

#funnel plot
Sedation.Funnel.Imp <- gg.funnel(es=Lab.Sedation.ES.Imp$es, es.var=Lab.Sedation.ES.Imp$var, 
                            mean.effect=rma.Sedation.Imp$ES.mean, se.effect=rma.Sedation.Imp$ES.SEM,
                            title="Lab Outcomes - Alcohol Sedation\nModerate Imputation", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                            lab=factor(Lab.Sedation.ES.Imp$Med), labsTitle="Medication")
Sedation.Funnel.Imp
#ggsave(Sedation.Funnel.Imp, filename="Sedation.Funnel.Imp.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Sedation.Imp$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = 0.5789, p = 0.5626

#Plot meta-analyzed effect sizes
rma.Sedation.Imp$ES.est$Med <- factor(rma.Sedation.Imp$ES.est$Med)
Sedation.ES.Imp.Plot <- ggplot(rma.Sedation.Imp$ES.est, aes(x=Se.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = Se.metaES - Se.metaES.se, xmax = Se.metaES + Se.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Sedation Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Sedation.Imp$ES.est$Med))) +
  ggtitle("Sedation Effect Sizes\nModerate Imputation") +
  SpTheme()
Sedation.ES.Imp.Plot
#ggsave(Sedation.ES.Imp.Plot, filename="Sedation.ES.Imp.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES.Imp <- full_join(Full.ES.Imp, rma.Sedation.Imp$ES.est, by="Med")



#NEGMOOD OUTCOME - Moderate Approach----

#impute no stat effect sizes to p=0.5
Lab.NegMood$ES.Imp <- NA
Lab.NegMood$ES.Impvar <- NA

for(i in 1:dim(Lab.NegMood)[1]){
  if(Lab.NegMood$NoStat[i]==1){
    Lab.NegMood$ES.Imp[i] <- as.double(pes(p=0.5, n.1=Lab.NegMood$N[i]/2, n.2=Lab.NegMood$N[i]/2, dig=4, verbose=F)[c("g")])
    Lab.NegMood$ES.Impvar[i] <- as.double(pes(p=0.5, n.1=Lab.NegMood$N[i]/2, n.2=Lab.NegMood$N[i]/2, dig=4, verbose=F)[c("var.g")])
  }
  else{
    Lab.NegMood$ES.Imp[i] <- Lab.NegMood$ES[i]
    Lab.NegMood$ES.Impvar[i] <- Lab.NegMood$ESvar[i]
  }
}

#checks
length(Lab.NegMood$ES.Imp) #46 effect sizes
length(table(Lab.NegMood$Med)) #12 medications with craving outcomes
table(Lab.NegMood$Sample) #which samples gave data
#Number of samples with craving data
length(levels(Lab.NegMood$Sample)) # 21 samples with craving outcomes

#Aggregate craving effect sizes
Lab.NegMood.ES.Imp <- agg(data=Lab.NegMood, id=Sample, es=ES.Imp, var=ES.Impvar,  method = "BHHR", cor=.6)
names(Lab.NegMood.ES.Imp)[names(Lab.NegMood.ES.Imp)=="id"] <- "Sample"
dim(Lab.NegMood.ES.Imp) #21 effect sizes

#merge aggregated effect sizes with 
dim(Lab.NegMood.ES.Imp)
dim(Lab.Samples)
Lab.NegMood.ES.Imp <- inner_join(Lab.NegMood.ES.Imp, Lab.Samples, by="Sample")
dim(Lab.NegMood.ES.Imp) #21 effect sizes across samples

#reset levels of Med and Sample
Lab.NegMood.ES.Imp$Med <- factor(Lab.NegMood.ES.Imp$Med)
Lab.NegMood.ES.Imp$Sample <- factor(Lab.NegMood.ES.Imp$Sample)

#count number of NegMood outcomes that were aggregated for each aggregated ES.Imp
#count outcomes
for (i in 1:dim(Lab.NegMood.ES.Imp)[1])
{
  Lab.NegMood.ES.Imp$Outcomes[i] <- nrow(subset(Lab.NegMood, Sample==Lab.NegMood.ES.Imp$Sample[i]))
}


#NegMood - RMA Analyses
rma.NegMood.Imp<- rma.recenterMed.Lab(Lab.NegMood.ES.Imp, abr="NM.")

#NegMood - Forest Plot
#Saving Size 10x10
forest.rma(rma.NegMood.Imp$rma.uncent,
           slab = paste(Lab.NegMood.ES.Imp$Author, Lab.NegMood.ES.Imp$Year,sep=", "),
           ilab = cbind(as.character(Lab.NegMood.ES.Imp$Med), Lab.NegMood.ES.Imp$MaxDose, round(Lab.NegMood.ES.Imp$DpM, 1), round(Lab.NegMood.ES.Imp$MaxAlcDose, 3)),
           ilab.xpos = c(-2.2, -1.2, 2, 2.7),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.NegMood.ES.Imp$Med),
           xlab="Hedge's G")
text(-4, 23, "Author(s) and Year", pos = 4, cex=1)
text(-2.2, 23, "Medication", pos = 4, cex=1)
text(-1.2, 23, "Dose", pos = 4, cex=1)
text(2, 23, "DpM", pos = 4, cex=1)
text(2.7, 23, "BrAC", pos = 4, cex=1)
text(3.3, 23, "Hedge's G [95% CI]", pos = 4, cex=1)
text(0, 24.5, "Alcohol Negative Mood - Moderate Imputation", cex=1.1)

#funnel plot
NegMood.Funnel.Imp <- gg.funnel(es=Lab.NegMood.ES.Imp$es, es.var=Lab.NegMood.ES.Imp$var, 
                             mean.effect=rma.NegMood.Imp$ES.mean, se.effect=rma.NegMood.Imp$ES.SEM,
                             title="Lab Outcomes - Alcohol NegMood\nModerate Imputation", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                             lab=factor(Lab.NegMood.ES.Imp$Med), labsTitle="Medication")
NegMood.Funnel.Imp <- NegMood.Funnel.Imp + annotate("rect", xmin = rma.NegMood.Imp$ES.mean+1.96*0.5, xmax = 1.3, ymin = 0, ymax = 0.5, alpha = .1, fill="black")
NegMood.Funnel.Imp
#ggsave(NegMood.Funnel.Imp, filename="NegMood.Funnel.Imp.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.NegMood.Imp$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = 1.0283, p = 0.3038

#Plot meta-analyzed effect sizes
rma.NegMood.Imp$ES.est$Med <- factor(rma.NegMood.Imp$ES.est$Med)
NegMood.ES.Imp.Plot <- ggplot(rma.NegMood.Imp$ES.est, aes(x=NM.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = NM.metaES - NM.metaES.se, xmax = NM.metaES + NM.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("NegMood Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.NegMood.Imp$ES.est$Med))) +
  ggtitle("NegMood Effect Sizes\nModerate Imputation") +
  SpTheme()
NegMood.ES.Imp.Plot
#ggsave(NegMood.ES.Imp.Plot, filename="NegMood.ES.Imp.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES.Imp <- full_join(Full.ES.Imp, rma.NegMood.Imp$ES.est, by="Med")


#HEAVY DRINKING OUTCOME - Moderate Approach----

#impute no stat effect sizes to p=0.5
RCT.Heavy$ES.Imp <- NA
RCT.Heavy$ES.Impvar <- NA

for(i in 1:dim(RCT.Heavy)[1]){
  if(RCT.Heavy$NoStat[i]==1){
    RCT.Heavy$ES.Imp[i] <- -as.double(pes(p=0.5, n.1=RCT.Heavy$N[i]/2, n.2=RCT.Heavy$N[i]/2, dig=4, verbose=F)[c("g")])
    RCT.Heavy$ES.Impvar[i] <- as.double(pes(p=0.5, n.1=RCT.Heavy$N[i]/2, n.2=RCT.Heavy$N[i]/2, dig=4, verbose=F)[c("var.g")])
  }
  else{
    RCT.Heavy$ES.Imp[i] <- RCT.Heavy$ES[i]
    RCT.Heavy$ES.Impvar[i] <- RCT.Heavy$ESvar[i]
  }
}

#checks
dim(RCT.Heavy) #406 effect sizes
table(RCT.Heavy$Med)
length(table(RCT.Heavy$Med)) #19 medications with Heavy Drinking outcomes
table(RCT.Heavy$Sample) #which samples gave data
#Number of samples with NegMood data
length(levels(RCT.Heavy$Sample)) # 132 samples with Heavy Drinking outcomes

#Aggregate Heavy effect sizes
RCT.Heavy.ES.Imp <- agg(data=RCT.Heavy, id=Sample, es=ES.Imp, var=ES.Impvar,  method = "BHHR", cor=.6)
names(RCT.Heavy.ES.Imp)[names(RCT.Heavy.ES.Imp)=="id"] <- "Sample"
dim(RCT.Heavy.ES.Imp) #132 effect sizes

#merge aggregated effect sizes with 
dim(RCT.Heavy.ES.Imp)
dim(Lab.Samples)
RCT.Heavy.ES.Imp <- inner_join(RCT.Heavy.ES.Imp, RCT.Samples, by="Sample")
dim(RCT.Heavy.ES.Imp) #132 effect sizes across samples

#reset levels of Med and Sample
RCT.Heavy.ES.Imp$Med <- factor(RCT.Heavy.ES.Imp$Med)
RCT.Heavy.ES.Imp$Sample <- factor(RCT.Heavy.ES.Imp$Sample)

#count number of Heavy outcomes that were aggregated for each aggregated ES.Imp
#count outcomes
for (i in 1:dim(RCT.Heavy.ES.Imp)[1])
{
  RCT.Heavy.ES.Imp$Outcomes[i] <- nrow(subset(RCT.Heavy, Sample==RCT.Heavy.ES.Imp$Sample[i]))
}

#Heavy - RMA Analyses
rma.Heavy.Imp<- rma.recenterMed.RCT(RCT.Heavy.ES.Imp, abr="He.")

#Heavy - Forest Plot
#Saving Size 10x15
forest.rma(rma.Heavy.Imp$rma.uncent,
           slab = paste(RCT.Heavy.ES.Imp$Author, RCT.Heavy.ES.Imp$Year,sep=", "),
           ilab = cbind(as.character(RCT.Heavy.ES.Imp$Med), RCT.Heavy.ES.Imp$MaxDose, round(RCT.Heavy.ES.Imp$TrxDur, 1)),
           ilab.xpos = c(-4, -3, 2.5),
           ilab.pos=c(4,4,4),
           order=order(RCT.Heavy.ES.Imp$Med),
           xlab="Hedge's G", 
           cex=0.5)
text(-5.5, 135, "Author(s) and Year", pos = 4, cex=0.6)
text(-4, 135, "Medication", pos = 4, cex=0.6)
text(-3, 135, "Dose", pos = 4, cex=0.6)
text(2.5, 135, "Treatment Duration", pos = 4, cex=0.6)
text(4, 135, "Hedge's G [95% CI]", pos = 4, cex=0.6)
text(0, 137, "Heavy Drinking - Moderate Imputation", cex=.8)

#funnel plot
Heavy.Funnel.Imp <- gg.funnel(es=RCT.Heavy.ES.Imp$es, es.var=RCT.Heavy.ES.Imp$var, 
                          mean.effect=rma.Heavy.Imp$ES.mean, se.effect=rma.Heavy.Imp$ES.SEM,
                          title="RCT Outcomes - Heavy Drinking\nModerate Imputation", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                          lab=factor(RCT.Heavy.ES.Imp$Med), labsTitle="Medication")
Heavy.Funnel.Imp
#ggsave(Heavy.Funnel.Imp, filename="Heavy.Funnel.Imp.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Heavy.Imp$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = -3.6156, p = 0.0003   #VERY SIGNIFICANT ASSYMETRY

#Plot meta-analyzed effect sizes
rma.Heavy.Imp$ES.Imp.est$Med <- factor(rma.Heavy.Imp$ES.est$Med)
Heavy.ES.Imp.Plot <- ggplot(rma.Heavy.Imp$ES.est, aes(x=He.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = He.metaES - He.metaES.se, xmax = He.metaES + He.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Heavy Drinking Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Heavy.Imp$ES.est$Med))) +
  ggtitle("Heavy Drinking Effect Sizes\nModerate Imputation") +
  SpTheme()
Heavy.ES.Imp.Plot
#ggsave(Heavy.ES.Imp.Plot, filename="Heavy.ES.Imp.Plot.png", width = 6, height = 5, dpi = 400)


#Save Medication Values
Full.ES.Imp <- full_join(Full.ES.Imp, rma.Heavy.Imp$ES.est, by="Med")



#ABSTINENCE DRINKING OUTCOME - Moderate Approach----

#impute no stat effect sizes to p=0.5
RCT.Abstinence$ES.Imp <- NA
RCT.Abstinence$ES.Impvar <- NA

for(i in 1:dim(RCT.Abstinence)[1]){
  if(RCT.Abstinence$NoStat[i]==1){
    RCT.Abstinence$ES.Imp[i] <- as.double(pes(p=0.5, n.1=RCT.Abstinence$N[i]/2, n.2=RCT.Abstinence$N[i]/2, dig=4, verbose=F)[c("g")])
    RCT.Abstinence$ES.Impvar[i] <- as.double(pes(p=0.5, n.1=RCT.Abstinence$N[i]/2, n.2=RCT.Abstinence$N[i]/2, dig=4, verbose=F)[c("var.g")])
  }
  else{
    RCT.Abstinence$ES.Imp[i] <- RCT.Abstinence$ES[i]
    RCT.Abstinence$ES.Impvar[i] <- RCT.Abstinence$ESvar[i]
  }
}

#checks
dim(RCT.Abstinence) #257 effect sizes
table(RCT.Abstinence$Med)
length(table(RCT.Abstinence$Med)) #19 medications with Abstinence Drinking outcomes
table(RCT.Abstinence$Sample) #which samples gave data
#Number of samples with NegMood data
length(levels(RCT.Abstinence$Sample)) # 129 samples with Abstinence Drinking outcomes

#Aggregate Abstinence effect sizes
RCT.Abstinence.ES.Imp <- agg(data=RCT.Abstinence, id=Sample, es=ES.Imp, var=ES.Impvar,  method = "BHHR", cor=.6)
names(RCT.Abstinence.ES.Imp)[names(RCT.Abstinence.ES.Imp)=="id"] <- "Sample"
dim(RCT.Abstinence.ES.Imp) #129 effect sizes

#merge aggregated effect sizes with 
dim(RCT.Abstinence.ES.Imp)
dim(Lab.Samples)
RCT.Abstinence.ES.Imp <- inner_join(RCT.Abstinence.ES.Imp, RCT.Samples, by="Sample")
dim(RCT.Abstinence.ES.Imp) #129 effect sizes across samples

#reset levels of Med and Sample
RCT.Abstinence.ES.Imp$Med <- factor(RCT.Abstinence.ES.Imp$Med)
RCT.Abstinence.ES.Imp$Sample <- factor(RCT.Abstinence.ES.Imp$Sample)

#count number of Abstinence outcomes that were aggregated for each aggregated ES.Imp
#count outcomes
for (i in 1:dim(RCT.Abstinence.ES.Imp)[1])
{
  RCT.Abstinence.ES.Imp$Outcomes[i] <- nrow(subset(RCT.Abstinence, Sample==RCT.Abstinence.ES.Imp$Sample[i]))
}

#Abstinence - RMA Analyses
rma.Abstinence.Imp<- rma.recenterMed.RCT(RCT.Abstinence.ES.Imp, abr="Ab.")

#Abstinence - Forest Plot
#Saving Size 10x15
forest.rma(rma.Abstinence.Imp$rma.uncent,
           slab = paste(RCT.Abstinence.ES.Imp$Author, RCT.Abstinence.ES.Imp$Year,sep=", "),
           ilab = cbind(as.character(RCT.Abstinence.ES.Imp$Med), RCT.Abstinence.ES.Imp$MaxDose, round(RCT.Abstinence.ES.Imp$TrxDur, 1)),
           ilab.xpos = c(-4, -2.5, 5),
           ilab.pos=c(4,4,4),
           order=order(RCT.Abstinence.ES.Imp$Med),
           xlab="Hedge's G", 
           cex=0.5)
text(-7.5, 131, "Author(s) and Year", pos = 4, cex=0.6)
text(-4, 131, "Medication", pos = 4, cex=0.6)
text(-2.5, 131, "Dose", pos = 4, cex=0.6)
text(5, 131, "Treatment Duration", pos = 4, cex=0.6)
text(8, 131, "Hedge's G [95% CI]", pos = 4, cex=0.6)
text(0, 133, "Abstinence - Moderate Imputation", cex=.8)

#funnel plot
Abstinence.Funnel.Imp <- gg.funnel(es=RCT.Abstinence.ES.Imp$es, es.var=RCT.Abstinence.ES.Imp$var, 
                              mean.effect=rma.Abstinence.Imp$ES.mean, se.effect=rma.Abstinence.Imp$ES.SEM,
                              title="RCT Outcomes - Abstinence\nModerate Imputation", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                              lab=factor(RCT.Abstinence.ES.Imp$Med), labsTitle="Medication")
Abstinence.Funnel.Imp <- Abstinence.Funnel.Imp + annotate("rect", xmin = rma.Abstinence.Imp$ES.mean+1.96*0.6, xmax = 3.2, ymin = 0, ymax = 0.6, alpha = .1, fill="black")
Abstinence.Funnel.Imp
#ggsave(Abstinence.Funnel.Imp, filename="Abstinence.Funnel.Imp.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Abstinence.Imp$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#Regression Test for Funnel Plot Asymmetry
# 
# model:     mixed-effects meta-regression model
# predictor: standard error
# 
# test for funnel plot asymmetry: z = 2.6568, p = 0.0079   #SIGNIFICANT ASSYMETRY

#Plot meta-analyzed effect sizes
rma.Abstinence.Imp$ES.Imp.est$Med <- factor(rma.Abstinence.Imp$ES.est$Med)
Abstinence.ES.Imp.Plot <- ggplot(rma.Abstinence.Imp$ES.est, aes(x=Ab.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = Ab.metaES - Ab.metaES.se, xmax = Ab.metaES + Ab.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Abstinence Drinking Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Abstinence.Imp$ES.est$Med))) +
  ggtitle("Abstinence Drinking Effect Sizes\nModerate Imputation") +
  SpTheme()
Abstinence.ES.Imp.Plot
#ggsave(Abstinence.ES.Imp.Plot, filename="Abstinence.ES.Imp.Plot.png", width = 6, height = 5, dpi = 400)


#Save Medication Values
Full.ES.Imp <- full_join(Full.ES.Imp, rma.Abstinence.Imp$ES.est, by="Med")




#WILLIAMSON-YORK TRANSLATIONAL ANALYSIS - Moderate Approach----

#Lab Craving - Heavy Drinking----
Cr.He.ES.Imp <- na.exclude(Full.ES.Imp[,c("Med", "Cr.metaES", "Cr.metaES.se", "He.metaES", "He.metaES.se")])
dim(Cr.He.ES.Imp)
#13 medications for this analysis

Cr.He.WYbwls.Imp <- WYbwls(x=Cr.He.ES.Imp$Cr.metaES, xsd=Cr.He.ES.Imp$Cr.metaES.se,
                           y=Cr.He.ES.Imp$He.metaES, ysd=Cr.He.ES.Imp$He.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 10.456 	 190.728 
# Slope 	 32.6999 	 587.656 
# 
# r:   	 0.1366 
# r^2: 	 0.0187 
# p:   	 0.95662279027573 

#Modify plot with specific labels
Cr.He.WYbwls.Imp.plot <- Cr.He.WYbwls.Imp$plot + 
  ggtitle("Laboratory Craving and RCT Heavy Drinking\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Craving (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)", limits=c(-0.65, 0.35))
Cr.He.WYbwls.Imp.plot
#ggsave(Cr.He.WYbwls.Imp.plot, filename="Cr.He.WYbwls.Imp.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
Cr.He.ES.Imp$meanSD <- (Cr.He.ES.Imp$Cr.metaES.se + Cr.He.ES.Imp$He.metaES.se)/2
Cr.He.ES.Imp$Wsize <- 1/(Cr.He.ES.Imp$meanSD^2)
Data.Ellipse <- CompEllipse(x=Cr.He.ES.Imp$Cr.metaES, xsd=Cr.He.ES.Imp$Cr.metaES.se, y=Cr.He.ES.Imp$He.metaES, ysd=Cr.He.ES.Imp$He.metaES.se)
Cr.He.Ellipse.Imp.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=Cr.He.ES.Imp, aes(x=Cr.metaES, y=He.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=0.25, y=0.2, label=paste("italic(R)[WY]^2 == ",round(Cr.He.WYbwls.Imp$r2, 3)), parse=TRUE) + 
  annotate("text", x=0.25, y=0.15, label=paste("italic(p)[WY] == ",round(Cr.He.WYbwls.Imp$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Craving and RCT Heavy Drinking\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Craving (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") + 
  SpTheme()
Cr.He.Ellipse.Imp.plot
#ggsave(Cr.He.Ellipse.Imp.plot, filename="Cr.He.Ellipse.Imp.plot.png", width = 6, height = 5, dpi = 400)


#Lab Stimulation - Heavy Drinking----
St.He.ES.Imp <- na.exclude(Full.ES.Imp[,c("Med", "St.metaES", "St.metaES.se", "He.metaES", "He.metaES.se")])
dim(St.He.ES.Imp)
#15 medications for this analysis

St.He.WYbwls.Imp <- WYbwls(x=St.He.ES.Imp$St.metaES, xsd=St.He.ES.Imp$St.metaES.se,
                           y=St.He.ES.Imp$He.metaES, ysd=St.He.ES.Imp$He.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#            Est  	   SE 
# Int   	 -0.198 	 0.13 
# Slope 	 1.6301 	 1.691 
# 
# r:   	 0.1156 
# r^2: 	 0.0134 
# p:   	 0.352726146774967  

#Modify plot with specific labels
St.He.WYbwls.Imp.plot <- St.He.WYbwls.Imp$plot + 
  ggtitle("Laboratory Stimulation and RCT Heavy Drinking\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Stimulation (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)", limits=c(-0.65, 0.35))
St.He.WYbwls.Imp.plot
#ggsave(St.He.WYbwls.Imp.plot, filename="St.He.WYbwls.Imp.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
St.He.ES.Imp$meanSD <- (St.He.ES.Imp$St.metaES.se + St.He.ES.Imp$He.metaES.se)/2
St.He.ES.Imp$Wsize <- 1/(St.He.ES.Imp$meanSD^2)
Data.Ellipse <- CompEllipse(x=St.He.ES.Imp$St.metaES, xsd=St.He.ES.Imp$St.metaES.se, y=St.He.ES.Imp$He.metaES, ysd=St.He.ES.Imp$He.metaES.se)
St.He.Ellipse.Imp.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=St.He.ES.Imp, aes(x=St.metaES, y=He.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=1, y=0.2, label=paste("italic(R)[WY]^2 == ",round(St.He.WYbwls.Imp$r2, 3)), parse=TRUE) + 
  annotate("text", x=1, y=0.15, label=paste("italic(p)[WY] == ",round(St.He.WYbwls.Imp$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Stimulation and RCT Heavy Drinking\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Stimulation (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") + 
  SpTheme()
St.He.Ellipse.Imp.plot
#ggsave(St.He.Ellipse.Imp.plot, filename="St.He.Ellipse.Imp.plot.png", width = 6, height = 5, dpi = 400)


#Lab Sedation - Heavy Drinking----
Se.He.ES.Imp <- na.exclude(Full.ES.Imp[,c("Med", "Se.metaES", "Se.metaES.se", "He.metaES", "He.metaES.se")])
dim(Se.He.ES.Imp)
#15 medications for this analysis

Se.He.WYbwls.Imp <- WYbwls(x=Se.He.ES.Imp$Se.metaES, xsd=Se.He.ES.Imp$Se.metaES.se,
                           y=Se.He.ES.Imp$He.metaES, ysd=Se.He.ES.Imp$He.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 0.256 	 0.91 
# Slope 	 -3.8446 	 7.969 
# 
# r:   	 -0.3151 
# r^2: 	 0.0993 
# p:   	 0.637522720399796 

#Modify plot with specific labels
Se.He.WYbwls.Imp.plot <- Se.He.WYbwls.Imp$plot + 
  ggtitle("Laboratory Sedation and RCT Heavy Drinking\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Sedation (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)", limits=c(-0.65, 0.35))
Se.He.WYbwls.Imp.plot
#ggsave(Se.He.WYbwls.Imp.plot, filename="Se.He.WYbwls.Imp.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
Se.He.ES.Imp$meanSD <- (Se.He.ES.Imp$Se.metaES.se + Se.He.ES.Imp$He.metaES.se)/2
Se.He.ES.Imp$Wsize <- 1/(Se.He.ES.Imp$meanSD^2)
Data.Ellipse <- CompEllipse(x=Se.He.ES.Imp$Se.metaES, xsd=Se.He.ES.Imp$Se.metaES.se, y=Se.He.ES.Imp$He.metaES, ysd=Se.He.ES.Imp$He.metaES.se)
Se.He.Ellipse.Imp.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=Se.He.ES.Imp, aes(x=Se.metaES, y=He.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=0.7, y=0.2, label=paste("italic(R)[WY]^2 == ",round(Se.He.WYbwls.Imp$r2, 3)), parse=TRUE) + 
  annotate("text", x=0.7, y=0.15, label=paste("italic(p)[WY] == ",round(Se.He.WYbwls.Imp$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Sedation and RCT Heavy Drinking\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Sedation (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") + 
  SpTheme()
Se.He.Ellipse.Imp.plot
#ggsave(Se.He.Ellipse.Imp.plot, filename="Se.He.Ellipse.Imp.plot.png", width = 6, height = 5, dpi = 400)


#Lab NegMood - Heavy Drinking----
NM.He.ES.Imp <- na.exclude(Full.ES.Imp[,c("Med", "NM.metaES", "NM.metaES.se", "He.metaES", "He.metaES.se")])
dim(NM.He.ES.Imp)
#15 medications for this analysis

NM.He.WYbwls.Imp <- WYbwls(x=NM.He.ES.Imp$NM.metaES, xsd=NM.He.ES.Imp$NM.metaES.se,
                           y=NM.He.ES.Imp$He.metaES, ysd=NM.He.ES.Imp$He.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 0.71 	 5.449 
# Slope 	 -4.7678 	 29.264 
# 
# r:   	 -0.3621 
# r^2: 	 0.1311 
# p:   	 0.875929733929411 

#Modify plot with specific labels
NM.He.WYbwls.Imp.plot <- NM.He.WYbwls.Imp$plot + 
  ggtitle("Laboratory Negative Mood and RCT Heavy Drinking\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Negative Mood (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)", limits=c(-0.65, 0.35))
NM.He.WYbwls.Imp.plot
#ggsave(NM.He.WYbwls.Imp.plot, filename="NM.He.WYbwls.Imp.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
NM.He.ES.Imp$meanSD <- (NM.He.ES.Imp$NM.metaES.se + NM.He.ES.Imp$He.metaES.se)/2
NM.He.ES.Imp$Wsize <- 1/(NM.He.ES.Imp$meanSD^2)
Data.Ellipse <- CompEllipse(x=NM.He.ES.Imp$NM.metaES, xsd=NM.He.ES.Imp$NM.metaES.se, y=NM.He.ES.Imp$He.metaES, ysd=NM.He.ES.Imp$He.metaES.se)
NM.He.Ellipse.Imp.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=NM.He.ES.Imp, aes(x=NM.metaES, y=He.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=-0.4, y=-0.3, label=paste("italic(R)[WY]^2 == ",round(NM.He.WYbwls.Imp$r2, 3)), parse=TRUE) + 
  annotate("text", x=-0.4, y=-0.35, label=paste("italic(p)[WY] == ",round(NM.He.WYbwls.Imp$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Negative Mood and RCT Heavy Drinking\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Negative Mood (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") + 
  SpTheme()
NM.He.Ellipse.Imp.plot
#ggsave(NM.He.Ellipse.Imp.plot, filename="NM.He.Ellipse.Imp.plot.png", width = 6, height = 5, dpi = 400)




#Lab Craving - Abstinence----
Cr.Ab.ES.Imp <- na.exclude(Full.ES.Imp[,c("Med", "Cr.metaES", "Cr.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(Cr.Ab.ES.Imp)
#13 medications for this analysis

Cr.Ab.WYbwls.Imp <- WYbwls(x=Cr.Ab.ES.Imp$Cr.metaES, xsd=Cr.Ab.ES.Imp$Cr.metaES.se,
                           y=Cr.Ab.ES.Imp$Ab.metaES, ysd=Cr.Ab.ES.Imp$Ab.metaES.se,
                           print=T, plot=T, tol=1e-08)

#  Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 0.527 	 0.385 
# Slope 	 1.3018 	 1.206 
# 
# r:   	 0.2171 
# r^2: 	 0.0471 
# p:   	 0.303509879010203 

#Modify plot with specific labels
Cr.Ab.WYbwls.Imp.plot <- Cr.Ab.WYbwls.Imp$plot + 
  ggtitle("Laboratory Craving and RCT Abstinence\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Craving (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
Cr.Ab.WYbwls.Imp.plot
#ggsave(Cr.Ab.WYbwls.Imp.plot, filename="Cr.Ab.WYbwls.Imp.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
Cr.Ab.ES.Imp$meanSD <- (Cr.Ab.ES.Imp$Cr.metaES.se + Cr.Ab.ES.Imp$Ab.metaES.se)/2
Cr.Ab.ES.Imp$Wsize <- 1/(Cr.Ab.ES.Imp$meanSD^2)
Data.Ellipse <- CompEllipse(x=Cr.Ab.ES.Imp$Cr.metaES, xsd=Cr.Ab.ES.Imp$Cr.metaES.se, y=Cr.Ab.ES.Imp$Ab.metaES, ysd=Cr.Ab.ES.Imp$Ab.metaES.se)
Cr.Ab.Ellipse.Imp.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=Cr.Ab.ES.Imp, aes(x=Cr.metaES, y=Ab.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=0.25, y=-0.3, label=paste("italic(R)[WY]^2 == ",round(Cr.Ab.WYbwls.Imp$r2, 3)), parse=TRUE) + 
  annotate("text", x=0.25, y=-0.35, label=paste("italic(p)[WY] == ",round(Cr.Ab.WYbwls.Imp$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Craving and RCT Abstinence\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Craving (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)") + 
  SpTheme()
Cr.Ab.Ellipse.Imp.plot
#ggsave(Cr.Ab.Ellipse.Imp.plot, filename="Cr.Ab.Ellipse.Imp.plot.png", width = 6, height = 5, dpi = 400)


#Lab Stimulation - Abstinence----
St.Ab.ES.Imp <- na.exclude(Full.ES.Imp[,c("Med", "St.metaES", "St.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(St.Ab.ES.Imp)
#15 medications for this analysis

St.Ab.WYbwls.Imp <- WYbwls(x=St.Ab.ES.Imp$St.metaES, xsd=St.Ab.ES.Imp$St.metaES.se,
                           y=St.Ab.ES.Imp$Ab.metaES, ysd=St.Ab.ES.Imp$Ab.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 0.135 	 0.033 
# Slope 	 -0.1316 	 0.123 
# 
# r:   	 -0.0499 
# r^2: 	 0.0025 
# p:   	 0.303096503815493 

#Modify plot with specific labels
St.Ab.WYbwls.Imp.plot <- St.Ab.WYbwls.Imp$plot + 
  ggtitle("Laboratory Stimulation and RCT Abstinence\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Stimulation (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
St.Ab.WYbwls.Imp.plot
#ggsave(St.Ab.WYbwls.Imp.plot, filename="St.Ab.WYbwls.Imp.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
St.Ab.ES.Imp$meanSD <- (St.Ab.ES.Imp$St.metaES.se + St.Ab.ES.Imp$Ab.metaES.se)/2
St.Ab.ES.Imp$Wsize <- 1/(St.Ab.ES.Imp$meanSD^2)
Data.Ellipse <- CompEllipse(x=St.Ab.ES.Imp$St.metaES, xsd=St.Ab.ES.Imp$St.metaES.se, y=St.Ab.ES.Imp$Ab.metaES, ysd=St.Ab.ES.Imp$Ab.metaES.se)
St.Ab.Ellipse.Imp.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=St.Ab.ES.Imp, aes(x=St.metaES, y=Ab.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=1, y=-0.2, label=paste("italic(R)[WY]^2 == ",round(St.Ab.WYbwls.Imp$r2, 3)), parse=TRUE) + 
  annotate("text", x=1, y=-0.25, label=paste("italic(p)[WY] == ",round(St.Ab.WYbwls.Imp$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Stimulation and RCT Abstinence\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Stimulation (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)") + 
  SpTheme()
St.Ab.Ellipse.Imp.plot
#ggsave(St.Ab.Ellipse.Imp.plot, filename="St.Ab.Ellipse.Imp.plot.png", width = 6, height = 5, dpi = 400)


#Lab Sedation - Abstinence----
Se.Ab.ES.Imp <- na.exclude(Full.ES.Imp[,c("Med", "Se.metaES", "Se.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(Se.Ab.ES.Imp)
#15 medications for this analysis

Se.Ab.WYbwls.Imp <- WYbwls(x=Se.Ab.ES.Imp$Se.metaES, xsd=Se.Ab.ES.Imp$Se.metaES.se,
                           y=Se.Ab.ES.Imp$Ab.metaES, ysd=Se.Ab.ES.Imp$Ab.metaES.se,
                           print=T, plot=T, tol=1e-08)

#  Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 0.1 	 0.05 
# Slope 	 0.3636 	 0.223 
# 
# r:   	 0.2632 
# r^2: 	 0.0693 
# p:   	 0.127747132637022 

#Modify plot with specific labels
Se.Ab.WYbwls.Imp.plot <- Se.Ab.WYbwls.Imp$plot + 
  ggtitle("Laboratory Sedation and RCT Abstinence\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Sedation (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
Se.Ab.WYbwls.Imp.plot
#ggsave(Se.Ab.WYbwls.Imp.plot, filename="Se.Ab.WYbwls.Imp.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
Se.Ab.ES.Imp$meanSD <- (Se.Ab.ES.Imp$Se.metaES.se + Se.Ab.ES.Imp$Ab.metaES.se)/2
Se.Ab.ES.Imp$Wsize <- 1/(Se.Ab.ES.Imp$meanSD^2)
Data.Ellipse <- CompEllipse(x=Se.Ab.ES.Imp$Se.metaES, xsd=Se.Ab.ES.Imp$Se.metaES.se, y=Se.Ab.ES.Imp$Ab.metaES, ysd=Se.Ab.ES.Imp$Ab.metaES.se)
Se.Ab.Ellipse.Imp.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=Se.Ab.ES.Imp, aes(x=Se.metaES, y=Ab.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=0.7, y=-0.4, label=paste("italic(R)[WY]^2 == ",round(Se.Ab.WYbwls.Imp$r2, 3)), parse=TRUE) + 
  annotate("text", x=0.7, y=-0.45, label=paste("italic(p)[WY] == ",round(Se.Ab.WYbwls.Imp$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Sedation and RCT Abstinence\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Sedation (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)") + 
  SpTheme()
Se.Ab.Ellipse.Imp.plot
#ggsave(Se.Ab.Ellipse.Imp.plot, filename="Se.Ab.Ellipse.Imp.plot.png", width = 6, height = 5, dpi = 400)


#Lab NegMood - Abstinence----
NM.Ab.ES.Imp <- na.exclude(Full.ES.Imp[,c("Med", "NM.metaES", "NM.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(NM.Ab.ES.Imp)
#8 medications for this analysis

NM.Ab.WYbwls.Imp <- WYbwls(x=NM.Ab.ES.Imp$NM.metaES, xsd=NM.Ab.ES.Imp$NM.metaES.se,
                           y=NM.Ab.ES.Imp$Ab.metaES, ysd=NM.Ab.ES.Imp$Ab.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 -0.346 	 1.442 
# Slope 	 2.6122 	 7.549 
# 
# r:   	 0.5755 
# r^2: 	 0.3312 
# p:   	 0.741111496067219 

#Modify plot with specific labels
NM.Ab.WYbwls.Imp.plot <- NM.Ab.WYbwls.Imp$plot + 
  ggtitle("Laboratory Negative Mood and RCT Abstinence\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Negative Mood (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
NM.Ab.WYbwls.Imp.plot
#ggsave(NM.Ab.WYbwls.Imp.plot, filename="NM.Ab.WYbwls.Imp.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
NM.Ab.ES.Imp$meanSD <- (NM.Ab.ES.Imp$NM.metaES.se + NM.Ab.ES.Imp$Ab.metaES.se)/2
NM.Ab.ES.Imp$Wsize <- 1/(NM.Ab.ES.Imp$meanSD^2)
Data.Ellipse <- CompEllipse(x=NM.Ab.ES.Imp$NM.metaES, xsd=NM.Ab.ES.Imp$NM.metaES.se, y=NM.Ab.ES.Imp$Ab.metaES, ysd=NM.Ab.ES.Imp$Ab.metaES.se)
NM.Ab.Ellipse.Imp.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=NM.Ab.ES.Imp, aes(x=NM.metaES, y=Ab.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=-0.4, y=0.28, label=paste("italic(R)[WY]^2 == ",round(NM.Ab.WYbwls.Imp$r2, 3)), parse=TRUE) + 
  annotate("text", x=-0.4, y=0.23, label=paste("italic(p)[WY] == ",round(NM.Ab.WYbwls.Imp$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Negative Mood and RCT Abstinence\n- Moderate Imputation -") +
  scale_x_continuous("Laboratory Effects on Alcohol Negative Mood (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)") + 
  SpTheme()
NM.Ab.Ellipse.Imp.plot
#ggsave(NM.Ab.Ellipse.Imp.plot, filename="NM.Ab.Ellipse.Imp.plot.png", width = 6, height = 5, dpi = 400)


#LIBERAL APROACH (No Stat = dropped)----

#CRAVING OUTCOME - Liberal Approach----

#drop no stat effect sizes
Lab.Craving.Drop <- subset(Lab.Craving, NoStat==0)

#reset levels of Med and Sample
Lab.Craving.Drop$Med <- factor(Lab.Craving.Drop$Med)
Lab.Craving.Drop$Sample <- factor(Lab.Craving.Drop$Sample)

#checks
length(Lab.Craving.Drop$ES) #49 effect sizes
length(table(Lab.Craving.Drop$Med)) #14 medications with craving outcomes
table(Lab.Craving.Drop$Sample) #which samples gave data
#Number of samples with craving data
length(levels(Lab.Craving.Drop$Sample)) # 32 samples with craving outcomes

#Aggregate craving effect sizes
Lab.Craving.Drop.ES <- agg(data=Lab.Craving.Drop, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(Lab.Craving.Drop.ES)[names(Lab.Craving.Drop.ES)=="id"] <- "Sample"
dim(Lab.Craving.Drop.ES) #32 effect sizes

#merge aggregated effect sizes with 
dim(Lab.Craving.Drop.ES)
dim(Lab.Samples)
Lab.Craving.Drop.ES <- inner_join(Lab.Craving.Drop.ES, Lab.Samples, by="Sample")
dim(Lab.Craving.Drop.ES) #32 effect sizes across samples

#reset levels of Med and Sample
Lab.Craving.Drop.ES$Med <- factor(Lab.Craving.Drop.ES$Med)
Lab.Craving.Drop.ES$Sample <- factor(Lab.Craving.Drop.ES$Sample)

#count number of Craving outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(Lab.Craving.Drop.ES)[1])
{
  Lab.Craving.Drop.ES$Outcomes[i] <- nrow(subset(Lab.Craving.Drop, Sample==Lab.Craving.Drop.ES$Sample[i]))
}

#Craving - RMA Analyses
rma.Craving.Drop<- rma.recenterMed.Lab(Lab.Craving.Drop.ES, abr="Cr.")

#Craving - Forest Plot
#Saving Size 10x10
forest.rma(rma.Craving.Drop$rma.uncent,
           slab = paste(Lab.Craving.Drop.ES$Author, Lab.Craving.Drop.ES$Year,sep=", "),
           ilab = cbind(as.character(Lab.Craving.Drop.ES$Med), Lab.Craving.Drop.ES$MaxDose, round(Lab.Craving.Drop.ES$DpM, 1), round(Lab.Craving.Drop.ES$MaxAlcDose, 3)),
           ilab.xpos = c(-3.3, -2.2, 1.2, 1.8),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.Craving.Drop.ES$Med),
           xlab="Hedge's G")
text(-5.3, 34, "Author(s) and Year", pos = 4, cex=1)
text(-3.3, 34, "Medication", pos = 4, cex=1)
text(-2.2, 34, "Dose", pos = 4, cex=1)
text(1.2, 34, "DpM", pos = 4, cex=1)
text(1.8, 34, "BrAC", pos = 4, cex=1)
text(2.6, 34, "Hedge's G [95% CI]", pos = 4, cex=1)
text(0, 36, "Alcohol Craving - Liberal Dropping", cex=1.2)

#funnel plot
Craving.Funnel.Drop <- gg.funnel(es=Lab.Craving.Drop.ES$es, es.var=Lab.Craving.Drop.ES$var, 
                            mean.effect=rma.Craving.Drop$ES.mean, se.effect=rma.Craving.Drop$ES.SEM,
                            title="Lab Outcomes - Alcohol Craving\nLiberal Dropping", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                            lab=factor(Lab.Craving.Drop.ES$Med), labsTitle="Medication")
Craving.Funnel.Drop
#ggsave(Craving.Funnel.Drop, filename="Craving.Funnel.Drop.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Craving.Drop$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#test for funnel plot asymmetry: z = -0.7210, p = 0.4709

#Plot meta-analyzed effect sizes
rma.Craving.Drop$ES.est$Med <- factor(rma.Craving.Drop$ES.est$Med)
Craving.ES.Drop.Plot <- ggplot(rma.Craving.Drop$ES.est, aes(x=Cr.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = Cr.metaES - Cr.metaES.se, xmax = Cr.metaES + Cr.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Craving Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Craving.Drop$ES.est$Med))) +
  ggtitle("Craving Effect Sizes\nLiberal Dropping") +
  SpTheme()
Craving.ES.Drop.Plot
#ggsave(Craving.ES.Drop.Plot, filename="Craving.ES.Drop.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES.Drop <- rma.Craving.Drop$ES.est


#STIMULATION OUTCOME - Liberal Approach----

#drop no stat effect sizes
Lab.Stimulation.Drop <- subset(Lab.Stimulation, NoStat==0)

#reset levels of Med and Sample
Lab.Stimulation.Drop$Med <- factor(Lab.Stimulation.Drop$Med)
Lab.Stimulation.Drop$Sample <- factor(Lab.Stimulation.Drop$Sample)

#checks
length(Lab.Stimulation.Drop$ES) #75 effect sizes
length(table(Lab.Stimulation.Drop$Med)) #17 medications with Stimulation outcomes
table(Lab.Stimulation.Drop$Sample) #which samples gave data
#Number of samples with craving data
length(levels(Lab.Stimulation.Drop$Sample)) # 38 samples with Stimulation outcomes

#Aggregate craving effect sizes
Lab.Stimulation.Drop.ES <- agg(data=Lab.Stimulation.Drop, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(Lab.Stimulation.Drop.ES)[names(Lab.Stimulation.Drop.ES)=="id"] <- "Sample"
dim(Lab.Stimulation.Drop.ES) #38 effect sizes

#merge aggregated effect sizes with 
dim(Lab.Stimulation.Drop.ES)
dim(Lab.Samples)
Lab.Stimulation.Drop.ES <- inner_join(Lab.Stimulation.Drop.ES, Lab.Samples, by="Sample")
dim(Lab.Stimulation.Drop.ES) #38 effect sizes across samples

#reset levels of Med and Sample
Lab.Stimulation.Drop.ES$Med <- factor(Lab.Stimulation.Drop.ES$Med)
Lab.Stimulation.Drop.ES$Sample <- factor(Lab.Stimulation.Drop.ES$Sample)

#count number of Stimulation outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(Lab.Stimulation.Drop.ES)[1])
{
  Lab.Stimulation.Drop.ES$Outcomes[i] <- nrow(subset(Lab.Stimulation.Drop, Sample==Lab.Stimulation.Drop.ES$Sample[i]))
}

#Stimulation - RMA Analyses
rma.Stimulation.Drop<- rma.recenterMed.Lab(Lab.Stimulation.Drop.ES, abr="St.")

#Stimulation - Forest Plot
#Saving Size 10x10
forest.rma(rma.Stimulation.Drop$rma.uncent,
           slab = paste(Lab.Stimulation.Drop.ES$Author, Lab.Stimulation.Drop.ES$Year,sep=", "),
           ilab = cbind(as.character(Lab.Stimulation.Drop.ES$Med), Lab.Stimulation.Drop.ES$MaxDose, round(Lab.Stimulation.Drop.ES$DpM, 1), round(Lab.Stimulation.Drop.ES$MaxAlcDose, 3)),
           ilab.xpos = c(-4, -2.5, 3.5, 4.5),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.Stimulation.Drop.ES$Med),
           xlab="Hedge's G")
text(-7, 40, "Author(s) and Year", pos = 4, cex=1)
text(-4, 40, "Medication", pos = 4, cex=1)
text(-2.5, 40, "Dose", pos = 4, cex=1)
text(3.5, 40, "DpM", pos = 4, cex=1)
text(4.5, 40, "BrAC", pos = 4, cex=1)
text(6, 40, "Hedge's G [95% CI]", pos = 4, cex=1)
text(0, 42, "Alcohol Stimulation - Liberal Dropping", cex=1.2)

#funnel plot
Stimulation.Funnel.Drop <- gg.funnel(es=Lab.Stimulation.Drop.ES$es, es.var=Lab.Stimulation.Drop.ES$var, 
                                 mean.effect=rma.Stimulation.Drop$ES.mean, se.effect=rma.Stimulation.Drop$ES.SEM,
                                 title="Lab Outcomes - Alcohol Stimulation\nLiberal Dropping", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                                 lab=factor(Lab.Stimulation.Drop.ES$Med), labsTitle="Medication")
Stimulation.Funnel.Drop <- Stimulation.Funnel.Drop  + annotate("rect", xmin = rma.Stimulation.Drop$ES.mean+1.96*0.6, xmax = 1.9, ymin = 0, ymax = 0.6, alpha = .1, fill="black")
Stimulation.Funnel.Drop
#ggsave(Stimulation.Funnel.Drop, filename="Stimulation.Funnel.Drop.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Stimulation.Drop$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#test for funnel plot asymmetry: z = -1.1149, p = 0.2649

#Plot meta-analyzed effect sizes
rma.Stimulation.Drop$ES.est$Med <- factor(rma.Stimulation.Drop$ES.est$Med)
Stimulation.ES.Drop.Plot <- ggplot(rma.Stimulation.Drop$ES.est, aes(x=St.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = St.metaES - St.metaES.se, xmax = St.metaES + St.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Stimulation Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Stimulation.Drop$ES.est$Med))) +
  ggtitle("Stimulation Effect Sizes\nLiberal Dropping") +
  SpTheme()
Stimulation.ES.Drop.Plot
#ggsave(Stimulation.ES.Drop.Plot, filename="Stimulation.ES.Drop.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES.Drop <- full_join(Full.ES.Drop, rma.Stimulation.Drop$ES.est, by="Med")


#SEDATION OUTCOME - Liberal Approach----

#drop no stat effect sizes
Lab.Sedation.Drop <- subset(Lab.Sedation, NoStat==0)

#reset levels of Med and Sample
Lab.Sedation.Drop$Med <- factor(Lab.Sedation.Drop$Med)
Lab.Sedation.Drop$Sample <- factor(Lab.Sedation.Drop$Sample)

#checks
length(Lab.Sedation.Drop$ES) #95 effect sizes
length(table(Lab.Sedation.Drop$Med)) #18 medications with Sedation outcomes
table(Lab.Sedation.Drop$Sample) #which samples gave data
#Number of samples with craving data
length(levels(Lab.Sedation.Drop$Sample)) # 42 samples with Sedation outcomes

#Aggregate craving effect sizes
Lab.Sedation.Drop.ES <- agg(data=Lab.Sedation.Drop, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(Lab.Sedation.Drop.ES)[names(Lab.Sedation.Drop.ES)=="id"] <- "Sample"
dim(Lab.Sedation.Drop.ES) #42 effect sizes

#merge aggregated effect sizes with 
dim(Lab.Sedation.Drop.ES)
dim(Lab.Samples)
Lab.Sedation.Drop.ES <- inner_join(Lab.Sedation.Drop.ES, Lab.Samples, by="Sample")
dim(Lab.Sedation.Drop.ES) #42 effect sizes across samples

#reset levels of Med and Sample
Lab.Sedation.Drop.ES$Med <- factor(Lab.Sedation.Drop.ES$Med)
Lab.Sedation.Drop.ES$Sample <- factor(Lab.Sedation.Drop.ES$Sample)

#count number of Sedation outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(Lab.Sedation.Drop.ES)[1])
{
  Lab.Sedation.Drop.ES$Outcomes[i] <- nrow(subset(Lab.Sedation.Drop, Sample==Lab.Sedation.Drop.ES$Sample[i]))
}

#Sedation - RMA Analyses
rma.Sedation.Drop<- rma.recenterMed.Lab(Lab.Sedation.Drop.ES, abr="Se.")

#Sedation - Forest Plot
#Saving Size 10x10
forest.rma(rma.Sedation.Drop$rma.uncent,
           slab = paste(Lab.Sedation.Drop.ES$Author, Lab.Sedation.Drop.ES$Year,sep=", "),
           ilab = cbind(as.character(Lab.Sedation.Drop.ES$Med), Lab.Sedation.Drop.ES$MaxDose, round(Lab.Sedation.Drop.ES$DpM, 1), round(Lab.Sedation.Drop.ES$MaxAlcDose, 3)),
           ilab.xpos = c(-3.3, -2.2, 2.5, 3.1),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.Sedation.Drop.ES$Med),
           xlab="Hedge's G")
text(-5.7, 44, "Author(s) and Year", pos = 4, cex=.8)
text(-3.3, 44, "Medication", pos = 4, cex=.8)
text(-2.2, 44, "Dose", pos = 4, cex=.8)
text(2.5, 44, "DpM", pos = 4, cex=.8)
text(3.1, 44, "BrAC", pos = 4, cex=.8)
text(3.9, 44, "Hedge's G [95% CI]", pos = 4, cex=.8)
text(0, 46, "Alcohol Sedation - Liberal Dropping", cex=1)

#funnel plot
Sedation.Funnel.Drop <- gg.funnel(es=Lab.Sedation.Drop.ES$es, es.var=Lab.Sedation.Drop.ES$var, 
                                 mean.effect=rma.Sedation.Drop$ES.mean, se.effect=rma.Sedation.Drop$ES.SEM,
                                 title="Lab Outcomes - Alcohol Sedation\nLiberal Dropping", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                                 lab=factor(Lab.Sedation.Drop.ES$Med), labsTitle="Medication")
Sedation.Funnel.Drop
#ggsave(Sedation.Funnel.Drop, filename="Sedation.Funnel.Drop.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Sedation.Drop$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#test for funnel plot asymmetry: z = -0.4852, p = 0.6275

#Plot meta-analyzed effect sizes
rma.Sedation.Drop$ES.est$Med <- factor(rma.Sedation.Drop$ES.est$Med)
Sedation.ES.Drop.Plot <- ggplot(rma.Sedation.Drop$ES.est, aes(x=Se.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = Se.metaES - Se.metaES.se, xmax = Se.metaES + Se.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Sedation Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Sedation.Drop$ES.est$Med))) +
  ggtitle("Sedation Effect Sizes\nLiberal Dropping") +
  SpTheme()
Sedation.ES.Drop.Plot
#ggsave(Sedation.ES.Drop.Plot, filename="Sedation.ES.Drop.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES.Drop <- full_join(Full.ES.Drop, rma.Sedation.Drop$ES.est, by="Med")


#NEGMOOD OUTCOME - Liberal Approach----

#drop no stat effect sizes
Lab.NegMood.Drop <- subset(Lab.NegMood, NoStat==0)

#reset levels of Med and Sample
Lab.NegMood.Drop$Med <- factor(Lab.NegMood.Drop$Med)
Lab.NegMood.Drop$Sample <- factor(Lab.NegMood.Drop$Sample)

#checks
length(Lab.NegMood.Drop$ES) #21 effect sizes
length(table(Lab.NegMood.Drop$Med)) #6 medications with NegMood outcomes
table(Lab.NegMood.Drop$Sample) #which samples gave data
#Number of samples with craving data
length(levels(Lab.NegMood.Drop$Sample)) #12 samples with NegMood outcomes

#Aggregate craving effect sizes
Lab.NegMood.Drop.ES <- agg(data=Lab.NegMood.Drop, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(Lab.NegMood.Drop.ES)[names(Lab.NegMood.Drop.ES)=="id"] <- "Sample"
dim(Lab.NegMood.Drop.ES) #12 effect sizes

#merge aggregated effect sizes with 
dim(Lab.NegMood.Drop.ES)
dim(Lab.Samples)
Lab.NegMood.Drop.ES <- inner_join(Lab.NegMood.Drop.ES, Lab.Samples, by="Sample")
dim(Lab.NegMood.Drop.ES) #12 effect sizes across samples

#reset levels of Med and Sample
Lab.NegMood.Drop.ES$Med <- factor(Lab.NegMood.Drop.ES$Med)
Lab.NegMood.Drop.ES$Sample <- factor(Lab.NegMood.Drop.ES$Sample)

#count number of NegMood outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(Lab.NegMood.Drop.ES)[1])
{
  Lab.NegMood.Drop.ES$Outcomes[i] <- nrow(subset(Lab.NegMood.Drop, Sample==Lab.NegMood.Drop.ES$Sample[i]))
}

#NegMood - RMA Analyses
rma.NegMood.Drop<- rma.recenterMed.Lab(Lab.NegMood.Drop.ES, abr="NM.")


#NegMood - Forest Plot
#Saving Size 10x10
forest.rma(rma.NegMood.Drop$rma.uncent,
           slab = paste(Lab.NegMood.Drop.ES$Author, Lab.NegMood.Drop.ES$Year,sep=", "),
           ilab = cbind(as.character(Lab.NegMood.Drop.ES$Med), Lab.NegMood.Drop.ES$MaxDose, round(Lab.NegMood.Drop.ES$DpM, 1), round(Lab.NegMood.Drop.ES$MaxAlcDose, 3)),
           ilab.xpos = c(-2.5, -1.6, 2.1, 2.6),
           ilab.pos=c(4,4,4,4),
           order=order(Lab.NegMood.Drop.ES$Med),
           xlab="Hedge's G")
text(-4, 14, "Author(s) and Year", pos = 4, cex=1)
text(-2.5, 14, "Medication", pos = 4, cex=1)
text(-1.6, 14, "Dose", pos = 4, cex=1)
text(2.1, 14, "DpM", pos = 4, cex=1)
text(2.6, 14, "BrAC", pos = 4, cex=1)
text(3.2, 14, "Hedge's G [95% CI]", pos = 4, cex=1)
text(0, 15, "Alcohol Negative Mood - Liberal Dropping", cex=1.2)

#funnel plot
NegMood.Funnel.Drop <- gg.funnel(es=Lab.NegMood.Drop.ES$es, es.var=Lab.NegMood.Drop.ES$var, 
                                 mean.effect=rma.NegMood.Drop$ES.mean, se.effect=rma.NegMood.Drop$ES.SEM,
                                 title="Lab Outcomes - Alcohol NegMood\nLiberal Dropping", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                                 lab=factor(Lab.NegMood.Drop.ES$Med), labsTitle="Medication")
NegMood.Funnel.Drop <- NegMood.Funnel.Drop  + annotate("rect", xmin = rma.NegMood.Drop$ES.mean+1.96*0.4, xmax = 1.2, ymin = 0, ymax = 0.4, alpha = .1, fill="black")
NegMood.Funnel.Drop
#ggsave(NegMood.Funnel.Drop, filename="NegMood.Funnel.Drop.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.NegMood.Drop$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#test for funnel plot asymmetry: z = 2.4136, p = 0.0158

#Plot meta-analyzed effect sizes
rma.NegMood.Drop$ES.est$Med <- factor(rma.NegMood.Drop$ES.est$Med)
NegMood.ES.Drop.Plot <- ggplot(rma.NegMood.Drop$ES.est, aes(x=NM.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = NM.metaES - NM.metaES.se, xmax = NM.metaES + NM.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("NegMood Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.NegMood.Drop$ES.est$Med))) +
  ggtitle("NegMood Effect Sizes\nLiberal Dropping") +
  SpTheme()
NegMood.ES.Drop.Plot
#ggsave(NegMood.ES.Drop.Plot, filename="NegMood.ES.Drop.Plot.png", width = 6, height = 5, dpi = 400)

#Save Medication Values
Full.ES.Drop <- full_join(Full.ES.Drop, rma.NegMood.Drop$ES.est, by="Med")


#HEAVY DRINKING OUTCOME - Liberal Approach----

#impute no stat effect sizes to p=0.5
RCT.Heavy.Drop <- subset(RCT.Heavy, NoStat==0)

#reset levels of Med and Sample
RCT.Heavy.Drop$Med <- factor(RCT.Heavy.Drop$Med)
RCT.Heavy.Drop$Sample <- factor(RCT.Heavy.Drop$Sample)

#checks
dim(RCT.Heavy.Drop) #182 effect sizes
table(RCT.Heavy.Drop$Med)
length(table(RCT.Heavy.Drop$Med)) #18 medications with Heavy Drinking outcomes
table(RCT.Heavy.Drop$Sample) #which samples gave data
#Number of samples with NegMood data
length(levels(RCT.Heavy.Drop$Sample)) # 104 samples with Heavy Drinking outcomes

#Aggregate Heavy effect sizes
RCT.Heavy.Drop.ES <- agg(data=RCT.Heavy.Drop, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(RCT.Heavy.Drop.ES)[names(RCT.Heavy.Drop.ES)=="id"] <- "Sample"
dim(RCT.Heavy.Drop.ES) #104 effect sizes

#merge aggregated effect sizes with 
dim(RCT.Heavy.Drop.ES)
dim(RCT.Samples)
RCT.Heavy.Drop.ES <- inner_join(RCT.Heavy.Drop.ES, RCT.Samples, by="Sample")
dim(RCT.Heavy.Drop.ES) #104 effect sizes across samples

#reset levels of Med and Sample
RCT.Heavy.Drop.ES$Med <- factor(RCT.Heavy.Drop.ES$Med)
RCT.Heavy.Drop.ES$Sample <- factor(RCT.Heavy.Drop.ES$Sample)

#count number of Heavy outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(RCT.Heavy.Drop.ES)[1])
{
  RCT.Heavy.Drop.ES$Outcomes[i] <- nrow(subset(RCT.Heavy.Drop, Sample==RCT.Heavy.Drop.ES$Sample[i]))
}

#Heavy - RMA Analyses
rma.Heavy.Drop<- rma.recenterMed.RCT(RCT.Heavy.Drop.ES, abr="He.")

#Heavy - Forest Plot
#Saving Size 10x15
forest.rma(rma.Heavy.Drop$rma.uncent,
           slab = paste(RCT.Heavy.Drop.ES$Author, RCT.Heavy.Drop.ES$Year,sep=", "),
           ilab = cbind(as.character(RCT.Heavy.Drop.ES$Med), RCT.Heavy.Drop.ES$MaxDose, round(RCT.Heavy.Drop.ES$TrxDur, 1)),
           ilab.xpos = c(-4, -3, 2.5),
           ilab.pos=c(4,4,4),
           order=order(RCT.Heavy.Drop.ES$Med),
           xlab="Hedge's G", 
           cex=0.5)
text(-6.5, 106, "Author(s) and Year", pos = 4, cex=0.6)
text(-4, 106, "Medication", pos = 4, cex=0.6)
text(-3, 106, "Dose", pos = 4, cex=0.6)
text(2.5, 106, "Treatment Duration", pos = 4, cex=0.6)
text(5.5, 106, "Hedge's G [95% CI]", pos = 4, cex=0.6)
text(0, 108, "Heavy Drinking - Liberal Dropping", cex=.8)

#funnel plot
Heavy.Funnel.Drop <- gg.funnel(es=RCT.Heavy.Drop.ES$es, es.var=RCT.Heavy.Drop.ES$var, 
                              mean.effect=rma.Heavy.Drop$ES.mean, se.effect=rma.Heavy.Drop$ES.SEM,
                              title="RCT Outcomes - Heavy Drinking\nLiberal Dropping", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                              lab=factor(RCT.Heavy.Drop.ES$Med), labsTitle="Medication")
Heavy.Funnel.Drop <- Heavy.Funnel.Drop + annotate("rect", xmin = rma.Heavy.Drop$ES.mean+1.96*0.6, xmax = 1.6, ymin = 0, ymax = 0.6, alpha = .1, fill="black")
Heavy.Funnel.Drop
#ggsave(Heavy.Funnel.Drop, filename="Heavy.Funnel.Drop.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Heavy.Drop$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#test for funnel plot asymmetry: z = -2.5376, p = 0.0112 #SIGNIFICANT ASSYMETRY

#Plot meta-analyzed effect sizes
rma.Heavy.Drop$ES.est$Med <- factor(rma.Heavy.Drop$ES.est$Med)
Heavy.ES.Drop.Plot <- ggplot(rma.Heavy.Drop$ES.est, aes(x=He.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = He.metaES - He.metaES.se, xmax = He.metaES + He.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Heavy Drinking Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Heavy.Drop$ES.est$Med))) +
  ggtitle("Heavy Drinking Effect Sizes\nLiberal Dropping") +
  SpTheme()
Heavy.ES.Drop.Plot
#ggsave(Heavy.ES.Drop.Plot, filename="Heavy.ES.Drop.Plot.png", width = 6, height = 5, dpi = 400)


#Save Medication Values
Full.ES.Drop <- full_join(Full.ES.Drop, rma.Heavy.Drop$ES.est, by="Med")


#ABSTINENCE DRINKING OUTCOME - Liberal Approach----

#impute no stat effect sizes to p=0.5
RCT.Abstinence.Drop <- subset(RCT.Abstinence, NoStat==0)

#reset levels of Med and Sample
RCT.Abstinence.Drop$Med <- factor(RCT.Abstinence.Drop$Med)
RCT.Abstinence.Drop$Sample <- factor(RCT.Abstinence.Drop$Sample)

#checks
dim(RCT.Abstinence.Drop) #154 effect sizes
table(RCT.Abstinence.Drop$Med)
length(table(RCT.Abstinence.Drop$Med)) #19 medications with Abstinence outcomes
table(RCT.Abstinence.Drop$Sample) #which samples gave data
#Number of samples with NegMood data
length(levels(RCT.Abstinence.Drop$Sample)) #111 samples with Abstinence outcomes

#Aggregate Abstinence effect sizes
RCT.Abstinence.Drop.ES <- agg(data=RCT.Abstinence.Drop, id=Sample, es=ES, var=ESvar,  method = "BHHR", cor=.6)
names(RCT.Abstinence.Drop.ES)[names(RCT.Abstinence.Drop.ES)=="id"] <- "Sample"
dim(RCT.Abstinence.Drop.ES) #111 effect sizes

#merge aggregated effect sizes with 
dim(RCT.Abstinence.Drop.ES)
dim(RCT.Samples)
RCT.Abstinence.Drop.ES <- inner_join(RCT.Abstinence.Drop.ES, RCT.Samples, by="Sample")
dim(RCT.Abstinence.Drop.ES) #111 effect sizes across samples

#reset levels of Med and Sample
RCT.Abstinence.Drop.ES$Med <- factor(RCT.Abstinence.Drop.ES$Med)
RCT.Abstinence.Drop.ES$Sample <- factor(RCT.Abstinence.Drop.ES$Sample)

#count number of Abstinence outcomes that were aggregated for each aggregated ES
#count outcomes
for (i in 1:dim(RCT.Abstinence.Drop.ES)[1])
{
  RCT.Abstinence.Drop.ES$Outcomes[i] <- nrow(subset(RCT.Abstinence.Drop, Sample==RCT.Abstinence.Drop.ES$Sample[i]))
}

#Abstinence - RMA Analyses
rma.Abstinence.Drop<- rma.recenterMed.RCT(RCT.Abstinence.Drop.ES, abr="Ab.")

#Abstinence - Forest Plot
#Saving Size 10x15
forest.rma(rma.Abstinence.Drop$rma.uncent,
           slab = paste(RCT.Abstinence.Drop.ES$Author, RCT.Abstinence.Drop.ES$Year,sep=", "),
           ilab = cbind(as.character(RCT.Abstinence.Drop.ES$Med), RCT.Abstinence.Drop.ES$MaxDose, round(RCT.Abstinence.Drop.ES$TrxDur, 1)),
           ilab.xpos = c(-6.5, -5, 5),
           ilab.pos=c(4,4,4),
           order=order(RCT.Abstinence.Drop.ES$Med),
           xlab="Hedge's G", 
           cex=0.5)
text(-9.5, 113, "Author(s) and Year", pos = 4, cex=0.6)
text(-6.5, 113, "Medication", pos = 4, cex=0.6)
text(-5, 113, "Dose", pos = 4, cex=0.6)
text(5, 113, "Treatment Duration", pos = 4, cex=0.6)
text(8.5, 113, "Hedge's G [95% CI]", pos = 4, cex=0.6)
text(0, 115, "Abstinence - Liberal Dropping", cex=.8)

#funnel plot
Abstinence.Funnel.Drop <- gg.funnel(es=RCT.Abstinence.Drop.ES$es, es.var=RCT.Abstinence.Drop.ES$var, 
                               mean.effect=rma.Abstinence.Drop$ES.mean, se.effect=rma.Abstinence.Drop$ES.SEM,
                               title="RCT Outcomes - Abstinence\nLiberal Dropping", x.lab="Effect Size (Hedge's G)", y.lab="Effect Size Std Error", 
                               lab=factor(RCT.Abstinence.Drop.ES$Med), labsTitle="Medication")
Abstinence.Funnel.Drop <- Abstinence.Funnel.Drop + annotate("rect", xmin = rma.Abstinence.Drop$ES.mean+1.96*0.6, xmax = 3.2, ymin = 0, ymax = 0.6, alpha = .1, fill="black") +
  annotate("rect", xmin = rma.Abstinence.Drop$ES.mean-1.96*0.6, xmax = -2, ymin = 0, ymax = 0.6, alpha = .1, fill="black")
Abstinence.Funnel.Drop
#ggsave(Abstinence.Funnel.Drop, filename="Abstinence.Funnel.Drop.png", width = 6, height = 5, dpi=400)

#test of funnel plot asymmetry
regtest(rma.Abstinence.Drop$rma.uncent, model="rma", predictor="sei", ret.fit=F)
#test for funnel plot asymmetry: z = -0.0490, p = 0.9609

#Plot meta-analyzed effect sizes
rma.Abstinence.Drop$ES.est$Med <- factor(rma.Abstinence.Drop$ES.est$Med)
Abstinence.ES.Drop.Plot <- ggplot(rma.Abstinence.Drop$ES.est, aes(x=Ab.metaES, y=Med)) +
  geom_vline(xintercept = 0, linetype='11') + 
  geom_errorbarh(aes(xmin = Ab.metaES - Ab.metaES.se, xmax = Ab.metaES + Ab.metaES.se), height = 0.2) + 
  geom_point(size=2) + 
  scale_x_continuous("Abstinence Effect Size (Hedge's g)") + 
  scale_y_discrete(name=element_blank(), limits=rev(levels(rma.Abstinence.Drop$ES.est$Med))) +
  ggtitle("Abstinence Effect Sizes\nLiberal Dropping") +
  SpTheme()
Abstinence.ES.Drop.Plot
#ggsave(Abstinence.ES.Drop.Plot, filename="Abstinence.ES.Drop.Plot.png", width = 6, height = 5, dpi = 400)


#Save Medication Values
Full.ES.Drop <- full_join(Full.ES.Drop, rma.Abstinence.Drop$ES.est, by="Med")



#WILLIAMSON-YORK TRANSLATIONAL ANALYSIS - Liberal Dropping----

#Lab Craving - Heavy Drinking----
Cr.He.ES.Drop <- na.exclude(Full.ES.Drop[,c("Med", "Cr.metaES", "Cr.metaES.se", "He.metaES", "He.metaES.se")])
dim(Cr.He.ES.Drop)
#11 medications for this analysis

Cr.He.WYbwlsDrop <- WYbwls(x=Cr.He.ES.Drop$Cr.metaES, xsd=Cr.He.ES.Drop$Cr.metaES.se,
                           y=Cr.He.ES.Drop$He.metaES, ysd=Cr.He.ES.Drop$He.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 -12.312 	 91.42 
# Slope 	 -32.028 	 239.517 
# 
# r:   	 -0.0067 
# r^2: 	 0 
# p:   	 0.896566832546067 

#Modify plot with specific labels
Cr.He.WYbwlsDrop.plot <- Cr.He.WYbwlsDrop$plot + 
  ggtitle("Laboratory Craving and RCT Heavy Drinking\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Craving (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)", limits=c(-1, 1.5))
Cr.He.WYbwlsDrop.plot
#ggsave(Cr.He.WYbwlsDrop.plot, filename="Cr.He.WYbwlsDrop.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
Cr.He.ES.Drop$meanSD <- (Cr.He.ES.Drop$Cr.metaES.se + Cr.He.ES.Drop$He.metaES.se)/2
Cr.He.ES.Drop$Wsize <- 1/(Cr.He.ES.Drop$meanSD^2)
Data.Ellipse <- CompEllipse(x=Cr.He.ES.Drop$Cr.metaES, xsd=Cr.He.ES.Drop$Cr.metaES.se, y=Cr.He.ES.Drop$He.metaES, ysd=Cr.He.ES.Drop$He.metaES.se)
Cr.He.EllipseDrop.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=Cr.He.ES.Drop, aes(x=Cr.metaES, y=He.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=0.25, y=1, label=paste("italic(R)[WY]^2 == ",round(Cr.He.WYbwlsDrop$r2, 3)), parse=TRUE) + 
  annotate("text", x=0.25, y=.9, label=paste("italic(p)[WY] == ",round(Cr.He.WYbwlsDrop$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Craving and RCT Heavy Drinking\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Craving (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") + 
  SpTheme()
Cr.He.EllipseDrop.plot
#ggsave(Cr.He.EllipseDrop.plot, filename="Cr.He.EllipseDrop.plot.png", width = 6, height = 5, dpi = 400)


#Lab Stimulation - Heavy Drinking----
St.He.ES.Drop <- na.exclude(Full.ES.Drop[,c("Med", "St.metaES", "St.metaES.se", "He.metaES", "He.metaES.se")])
dim(St.He.ES.Drop)
#12 medications for this analysis

St.He.WYbwlsDrop <- WYbwls(x=St.He.ES.Drop$St.metaES, xsd=St.He.ES.Drop$St.metaES.se,
                           y=St.He.ES.Drop$He.metaES, ysd=St.He.ES.Drop$He.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 -0.451 	 0.253 
# Slope 	 2.1876 	 1.621 
# 
# r:   	 0.1636 
# r^2: 	 0.0268 
# p:   	 0.206941761309258 

#Modify plot with specific labels
St.He.WYbwlsDrop.plot <- St.He.WYbwlsDrop$plot + 
  ggtitle("Laboratory Stimulation and RCT Heavy Drinking\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Stimulation (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)", limits=c(-1.1, 2.5)) + theme(legend.position=c(.25,.85))
St.He.WYbwlsDrop.plot
#ggsave(St.He.WYbwlsDrop.plot, filename="St.He.WYbwlsDrop.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
St.He.ES.Drop$meanSD <- (St.He.ES.Drop$St.metaES.se + St.He.ES.Drop$He.metaES.se)/2
St.He.ES.Drop$Wsize <- 1/(St.He.ES.Drop$meanSD^2)
Data.Ellipse <- CompEllipse(x=St.He.ES.Drop$St.metaES, xsd=St.He.ES.Drop$St.metaES.se, y=St.He.ES.Drop$He.metaES, ysd=St.He.ES.Drop$He.metaES.se)
St.He.EllipseDrop.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=St.He.ES.Drop, aes(x=St.metaES, y=He.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=1, y=1, label=paste("italic(R)[WY]^2 == ",round(St.He.WYbwlsDrop$r2, 3)), parse=TRUE) + 
  annotate("text", x=1, y=.9, label=paste("italic(p)[WY] == ",round(St.He.WYbwlsDrop$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Stimulation and RCT Heavy Drinking\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Stimulation (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") + 
  SpTheme()
St.He.EllipseDrop.plot
#ggsave(St.He.EllipseDrop.plot, filename="St.He.EllipseDrop.plot.png", width = 6, height = 5, dpi = 400)


#Lab Sedation - Heavy Drinking----
Se.He.ES.Drop <- na.exclude(Full.ES.Drop[,c("Med", "Se.metaES", "Se.metaES.se", "He.metaES", "He.metaES.se")])
dim(Se.He.ES.Drop)
#12 medications for this analysis

Se.He.WYbwlsDrop <- WYbwls(x=Se.He.ES.Drop$Se.metaES, xsd=Se.He.ES.Drop$Se.metaES.se,
                           y=Se.He.ES.Drop$He.metaES, ysd=Se.He.ES.Drop$He.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 0.249 	 0.833 
# Slope 	 -6.0245 	 9.055 
# 
# r:   	 -0.1896 
# r^2: 	 0.0359 
# p:   	 0.520893715067915 

#Modify plot with specific labels
Se.He.WYbwlsDrop.plot <- Se.He.WYbwlsDrop$plot + 
  ggtitle("Laboratory Sedation and RCT Heavy Drinking\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Sedation (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") + theme(legend.position=c(.25,.3))
Se.He.WYbwlsDrop.plot
#ggsave(Se.He.WYbwlsDrop.plot, filename="Se.He.WYbwlsDrop.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
Se.He.ES.Drop$meanSD <- (Se.He.ES.Drop$Se.metaES.se + Se.He.ES.Drop$He.metaES.se)/2
Se.He.ES.Drop$Wsize <- 1/(Se.He.ES.Drop$meanSD^2)
Data.Ellipse <- CompEllipse(x=Se.He.ES.Drop$Se.metaES, xsd=Se.He.ES.Drop$Se.metaES.se, y=Se.He.ES.Drop$He.metaES, ysd=Se.He.ES.Drop$He.metaES.se)
Se.He.EllipseDrop.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=Se.He.ES.Drop, aes(x=Se.metaES, y=He.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=.8, y=1, label=paste("italic(R)[WY]^2 == ",round(Se.He.WYbwlsDrop$r2, 3)), parse=TRUE) + 
  annotate("text", x=.8, y=.9, label=paste("italic(p)[WY] == ",round(Se.He.WYbwlsDrop$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Sedation and RCT Heavy Drinking\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Sedation (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)") + 
  SpTheme()
Se.He.EllipseDrop.plot
#ggsave(Se.He.EllipseDrop.plot, filename="Se.He.EllipseDrop.plot.png", width = 6, height = 5, dpi = 400)



#Lab NegMood - Heavy Drinking----
NM.He.ES.Drop <- na.exclude(Full.ES.Drop[,c("Med", "NM.metaES", "NM.metaES.se", "He.metaES", "He.metaES.se")])
dim(NM.He.ES.Drop)
#4 medications for this analysis

NM.He.WYbwlsDrop <- WYbwls(x=NM.He.ES.Drop$NM.metaES, xsd=NM.He.ES.Drop$NM.metaES.se,
                           y=NM.He.ES.Drop$He.metaES, ysd=NM.He.ES.Drop$He.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 -0.803 	 4.072 
# Slope 	 2.4912 	 19.663 
# 
# r:   	 0.3154 
# r^2: 	 0.0995 
# p:   	 0.910770197807318 

#Modify plot with specific labels
NM.He.WYbwlsDrop.plot <- NM.He.WYbwlsDrop$plot + 
  ggtitle("Laboratory Negative Mood and RCT Heavy Drinking\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol NegMood (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)")
NM.He.WYbwlsDrop.plot
#ggsave(NM.He.WYbwlsDrop.plot, filename="NM.He.WYbwlsDrop.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
NM.He.ES.Drop$meanSD <- (NM.He.ES.Drop$NM.metaES.se + NM.He.ES.Drop$He.metaES.se)/2
NM.He.ES.Drop$Wsize <- 1/(NM.He.ES.Drop$meanSD^2)
Data.Ellipse <- CompEllipse(x=NM.He.ES.Drop$NM.metaES, xsd=NM.He.ES.Drop$NM.metaES.se, y=NM.He.ES.Drop$He.metaES, ysd=NM.He.ES.Drop$He.metaES.se)
NM.He.EllipseDrop.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=NM.He.ES.Drop, aes(x=NM.metaES, y=He.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=0.7, y=-0.6, label=paste("italic(R)[WY]^2 == ",round(NM.He.WYbwlsDrop$r2, 3)), parse=TRUE) + 
  annotate("text", x=0.7, y=-0.65, label=paste("italic(p)[WY] == ",round(NM.He.WYbwlsDrop$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Negative Mood and RCT Heavy Drinking\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol NegMood (Hedge's G)") +
  scale_y_continuous("RCT Heavy Drinking Outcomes (Hedge's G)", limits = c(-0.7, 0.1)) + 
  SpTheme()
NM.He.EllipseDrop.plot
#ggsave(NM.He.EllipseDrop.plot, filename="NM.He.EllipseDrop.plot.png", width = 6, height = 5, dpi = 400)




#Lab Craving - Abstinence----
Cr.Ab.ES.Drop <- na.exclude(Full.ES.Drop[,c("Med", "Cr.metaES", "Cr.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(Cr.Ab.ES.Drop)
#11 medications for this analysis

Cr.Ab.WYbwlsDrop <- WYbwls(x=Cr.Ab.ES.Drop$Cr.metaES, xsd=Cr.Ab.ES.Drop$Cr.metaES.se,
                           y=Cr.Ab.ES.Drop$Ab.metaES, ysd=Cr.Ab.ES.Drop$Ab.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 1.954 	 2.815 
# Slope 	 5.1363 	 7.326 
# 
# r:   	 0.2345 
# r^2: 	 0.055 
# p:   	 0.500963589983834 

#Modify plot with specific labels
Cr.Ab.WYbwlsDrop.plot <- Cr.Ab.WYbwlsDrop$plot + 
  ggtitle("Laboratory Craving and RCT Abstinence\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Craving (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
Cr.Ab.WYbwlsDrop.plot
#ggsave(Cr.Ab.WYbwlsDrop.plot, filename="Cr.Ab.WYbwlsDrop.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
Cr.Ab.ES.Drop$meanSD <- (Cr.Ab.ES.Drop$Cr.metaES.se + Cr.Ab.ES.Drop$Ab.metaES.se)/2
Cr.Ab.ES.Drop$Wsize <- 1/(Cr.Ab.ES.Drop$meanSD^2)
Data.Ellipse <- CompEllipse(x=Cr.Ab.ES.Drop$Cr.metaES, xsd=Cr.Ab.ES.Drop$Cr.metaES.se, y=Cr.Ab.ES.Drop$Ab.metaES, ysd=Cr.Ab.ES.Drop$Ab.metaES.se)
Cr.Ab.EllipseDrop.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=Cr.Ab.ES.Drop, aes(x=Cr.metaES, y=Ab.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=0.25, y=1, label=paste("italic(R)[WY]^2 == ",round(Cr.Ab.WYbwlsDrop$r2, 3)), parse=TRUE) + 
  annotate("text", x=0.25, y=.9, label=paste("italic(p)[WY] == ",round(Cr.Ab.WYbwlsDrop$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Craving and RCT Abstinence\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Craving (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)") + 
  SpTheme()
Cr.Ab.EllipseDrop.plot
#ggsave(Cr.Ab.EllipseDrop.plot, filename="Cr.Ab.EllipseDrop.plot.png", width = 6, height = 5, dpi = 400)


#Lab Stimulation - Abstinence----
St.Ab.ES.Drop <- na.exclude(Full.ES.Drop[,c("Med", "St.metaES", "St.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(St.Ab.ES.Drop)
#12 medications for this analysis

St.Ab.WYbwlsDrop <- WYbwls(x=St.Ab.ES.Drop$St.metaES, xsd=St.Ab.ES.Drop$St.metaES.se,
                           y=St.Ab.ES.Drop$Ab.metaES, ysd=St.Ab.ES.Drop$Ab.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 0.165 	 0.092 
# Slope 	 -0.6643 	 0.323 
# 
# r:   	 -0.1564 
# r^2: 	 0.0245 
# p:   	 0.0669699255622393        #TREND LEVEL EFFECT GREATER REDUCTIONS IN STIMULATION ASSOCIATED WITH GREATER ABSTINENCE

#Modify plot with specific labels
St.Ab.WYbwlsDrop.plot <- St.Ab.WYbwlsDrop$plot + 
  ggtitle("Laboratory Stimulation and RCT Abstinence\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Stimulation (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
St.Ab.WYbwlsDrop.plot
#ggsave(St.Ab.WYbwlsDrop.plot, filename="St.Ab.WYbwlsDrop.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
St.Ab.ES.Drop$meanSD <- (St.Ab.ES.Drop$St.metaES.se + St.Ab.ES.Drop$Ab.metaES.se)/2
St.Ab.ES.Drop$Wsize <- 1/(St.Ab.ES.Drop$meanSD^2)
Data.Ellipse <- CompEllipse(x=St.Ab.ES.Drop$St.metaES, xsd=St.Ab.ES.Drop$St.metaES.se, y=St.Ab.ES.Drop$Ab.metaES, ysd=St.Ab.ES.Drop$Ab.metaES.se)
St.Ab.EllipseDrop.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=St.Ab.ES.Drop, aes(x=St.metaES, y=Ab.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=1.2, y=-0.8, label=paste("italic(R)[WY]^2 == ",round(St.Ab.WYbwlsDrop$r2, 3)), parse=TRUE) + 
  annotate("text", x=1.2, y=-0.9, label=paste("italic(p)[WY] == ",round(St.Ab.WYbwlsDrop$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Stimulation and RCT Abstinence\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Stimulation (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)") + 
  SpTheme()
St.Ab.EllipseDrop.plot
#ggsave(St.Ab.EllipseDrop.plot, filename="St.Ab.EllipseDrop.plot.png", width = 6, height = 5, dpi = 400)


#Lab Sedation - Abstinence----
Se.Ab.ES.Drop <- na.exclude(Full.ES.Drop[,c("Med", "Se.metaES", "Se.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(Se.Ab.ES.Drop)
#12 medications for this analysis

Se.Ab.WYbwlsDrop <- WYbwls(x=Se.Ab.ES.Drop$Se.metaES, xsd=Se.Ab.ES.Drop$Se.metaES.se,
                           y=Se.Ab.ES.Drop$Ab.metaES, ysd=Se.Ab.ES.Drop$Ab.metaES.se,
                           print=T, plot=T, tol=1e-08)

# Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 0.023 	 0.103 
# Slope 	 0.8639 	 0.454 
# 
# r:   	 0.3061 
# r^2: 	 0.0937 
# p:   	 0.086137515075164 

#Modify plot with specific labels
Se.Ab.WYbwlsDrop.plot <- Se.Ab.WYbwlsDrop$plot + 
  ggtitle("Laboratory Sedation and RCT Abstinence\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Sedation (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
Se.Ab.WYbwlsDrop.plot
#ggsave(Se.Ab.WYbwlsDrop.plot, filename="Se.Ab.WYbwlsDrop.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
Se.Ab.ES.Drop$meanSD <- (Se.Ab.ES.Drop$Se.metaES.se + Se.Ab.ES.Drop$Ab.metaES.se)/2
Se.Ab.ES.Drop$Wsize <- 1/(Se.Ab.ES.Drop$meanSD^2)
Data.Ellipse <- CompEllipse(x=Se.Ab.ES.Drop$Se.metaES, xsd=Se.Ab.ES.Drop$Se.metaES.se, y=Se.Ab.ES.Drop$Ab.metaES, ysd=Se.Ab.ES.Drop$Ab.metaES.se)
Se.Ab.EllipseDrop.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=Se.Ab.ES.Drop, aes(x=Se.metaES, y=Ab.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=.8, y=-0.6, label=paste("italic(R)[WY]^2 == ",round(Se.Ab.WYbwlsDrop$r2, 3)), parse=TRUE) + 
  annotate("text", x=.8, y=-0.7, label=paste("italic(p)[WY] == ",round(Se.Ab.WYbwlsDrop$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Sedation and RCT Abstinence\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol Sedation (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)") + 
  SpTheme()
Se.Ab.EllipseDrop.plot
#ggsave(Se.Ab.EllipseDrop.plot, filename="Se.Ab.EllipseDrop.plot.png", width = 6, height = 5, dpi = 400)



#Lab NegMood - Abstinence----
NM.Ab.ES.Drop <- na.exclude(Full.ES.Drop[,c("Med", "NM.metaES", "NM.metaES.se", "Ab.metaES", "Ab.metaES.se")])
dim(NM.Ab.ES.Drop)
#4 medications for this analysis

NM.Ab.WYbwlsDrop <- WYbwls(x=NM.Ab.ES.Drop$NM.metaES, xsd=NM.Ab.ES.Drop$NM.metaES.se,
                           y=NM.Ab.ES.Drop$Ab.metaES, ysd=NM.Ab.ES.Drop$Ab.metaES.se,
                           print=T, plot=T, tol=1e-08)

# # Williamson-York Algorithm for Bivariate Weighted Least Squared 
# 
# Coefficients: 
#   Est  	   SE 
# Int   	 -0.176 	 1.551 
# Slope 	 1.9188 	 7.195 
# 
# r:   	 0.3558 
# r^2: 	 0.1266 
# p:   	 0.814691624182437 

#Modify plot with specific labels
NM.Ab.WYbwlsDrop.plot <- NM.Ab.WYbwlsDrop$plot + 
  ggtitle("Laboratory Negative Mood and RCT Abstinence\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol NegMood (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)")
NM.Ab.WYbwlsDrop.plot
#ggsave(NM.Ab.WYbwlsDrop.plot, filename="NM.Ab.WYbwlsDrop.plot.png", width = 6, height = 5, dpi = 400)

#Plot without the regression lines
NM.Ab.ES.Drop$meanSD <- (NM.Ab.ES.Drop$NM.metaES.se + NM.Ab.ES.Drop$Ab.metaES.se)/2
NM.Ab.ES.Drop$Wsize <- 1/(NM.Ab.ES.Drop$meanSD^2)
Data.Ellipse <- CompEllipse(x=NM.Ab.ES.Drop$NM.metaES, xsd=NM.Ab.ES.Drop$NM.metaES.se, y=NM.Ab.ES.Drop$Ab.metaES, ysd=NM.Ab.ES.Drop$Ab.metaES.se)
NM.Ab.EllipseDrop.plot <- ggplot() +
  geom_hline(yintercept = 0, linetype='33') +
  geom_vline(xintercept = 0, linetype='33') +
  geom_polygon(data=Data.Ellipse,aes(x=xEll,y=yEll, group=obs), alpha=.15) +
  geom_point(data=NM.Ab.ES.Drop, aes(x=NM.metaES, y=Ab.metaES, size=Wsize), show.legend=F) +
  annotate("text", x=0.8, y=0.5, label=paste("italic(R)[WY]^2 == ",round(NM.Ab.WYbwlsDrop$r2, 3)), parse=TRUE) + 
  annotate("text", x=0.8, y=0.45, label=paste("italic(p)[WY] == ",round(NM.Ab.WYbwlsDrop$p, 3)), parse=TRUE) + 
  scale_size_continuous(range = c(2,7)) +
  ggtitle("Laboratory Negative Mood and RCT Abstinence\n- Liberal Dropping -") +
  scale_x_continuous("Laboratory Effects on Alcohol NegMood (Hedge's G)") +
  scale_y_continuous("RCT Abstinence Outcomes (Hedge's G)") + 
  SpTheme()
NM.Ab.EllipseDrop.plot
#ggsave(NM.Ab.EllipseDrop.plot, filename="NM.Ab.EllipseDrop.plot.png", width = 6, height = 5, dpi = 400)



