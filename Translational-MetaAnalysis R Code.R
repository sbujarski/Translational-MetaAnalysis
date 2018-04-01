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
library(xlsx) #package to import xls files directly



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


#IMPORT DATA----

Lab <- read.xlsx()

