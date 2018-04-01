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
library()


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
      #geom_polygon(data=CI95shade, aes(x=x, y=y), alpha=.3)+
      geom_vline(xintercept=mean.effect, size=1)+
      geom_vline(xintercept=0, linetype="dotted", size=1)+
      geom_point(size=3, aes(colour=lab))+
      geom_text(size=3, hjust = 0, nudge_x = 0.03, check_overlap = TRUE, aes(label=abbreviate(lab, minlength = 3)))+
      ggtitle(title)+
      xlab(x.lab)+
      ylab(y.lab)+
      theme_bw()+
      scale_y_reverse()+
      guides(colour=guide_legend(title=labsName))
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
      geom_polygon(data=CI95shade, aes(x=x, y=y), alpha=.3)+
      geom_vline(xintercept=mean.effect, size=1)+
      geom_vline(xintercept=0, linetype="dotted", size=1)+
      geom_point(size=3)+
      ggtitle(title)+
      xlab(x.lab)+
      ylab(y.lab)+
      theme_bw()+
      scale_y_reverse()
  }
}

#rma.recenterMed.Lab
#Function takes a pre-processed dataset with Med, es, var and covariate values for Laboratory data
#runs random effects meta-analysis with med, MaxDose.C (log modal centered medication dose), MaxAlcDose.C (0.06 centered), and DPM.C (gand mean centered) covariates
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
  cat("Random Effects Meta-Analysis result", "\n")
  rma.uncent <- rma(yi=es, vi=var, mods = ~ Med  + MaxDose.C + MaxAlcDose.C + logDpM.C, #meta regression outcome
                    data = data, method="REML")
  print(rma.uncent)
  
  #Get results averaging accross meds for funnel plot
  rma.noMed <- rma(yi=es, vi=var, mods = ~ MaxDose.C + MaxAlcDose.C + logDpM.C, data = data, method="REML")
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
    rma.result <- rma(yi=es, vi=var, mods = ~ Med  + MaxDose.C + MaxAlcDose.C + logDpM.C, #meta regression outcome
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

#CLEAN DATA----

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
ggsave(Lab.Years.Hist, filename="Lab.Years.Hist.png", width = 6, height = 5, dpi=300)

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
ggsave(Lab.N.Hist, filename="Lab.N.Hist.png", width = 6, height = 5, dpi=300)

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
ggsave(Lab.DpM.Hist, filename="Lab.DpM.Hist.png", width = 6, height = 5, dpi=300)
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
ggsave(Lab.MaxAlcDose.Hist, filename="Lab.MaxAlcDose.Hist.png", width = 6, height = 5, dpi=300)
#qq plot to look for normality of DpM
SpQQPlot(Lab.Samples$MaxAlcDose)
shapiro.test(Lab.Samples$MaxAlcDose) #p < 0.0001
SpQQPlot(log(Lab.Samples$MaxAlcDose))
shapiro.test(log(Lab.Samples$DpM)) #p = 0.06541


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

#Outcome domains
table(Lab.noSA.main$OutDomain)
table(Lab.noSA.main$OutDomain)/408
# Craving     NegMood    Sedation Stimulation 
#      72          46         171         119 
#     18%         11%         42%         29%
OutDomain.plot <- ggplot(Lab.noSA.main, aes(OutDomain)) + geom_bar(aes(fill=as.factor(NoStat)), width = 0.8) + 
  ggtitle("Number of Effect Sizes Per Domain") + 
  scale_x_discrete("Outcome Domain") +
  SpTheme(legend.position = "right")
OutDomain.plot



#IMPORT RCT DATA----
RCT <- read_excel("C:/Users/sbuja/Documents/Manuscripts for Publication/Quant Review Trans Addict Med/Meta Analysis/RCT Meta-Analysis Scored Data.xlsx", 
                  sheet = "Study Data", na="NA")


