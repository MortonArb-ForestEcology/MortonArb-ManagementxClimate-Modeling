#----------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick, Christy Rollinson
# Project: Mandifore -- Morton Arb Management Case Study
# Purpose: This script performs a cleaned up version of our Superimposed epoch analysis
# Inputs: Yearly ED2 output csv from script 2a_Yearly_ED2_sum.R
# Outputs: Figures and Tables
# Notes: This is a cleaner version of the anlysis scripts found in 6a_SEA_analysis.R
#----------------------------------------------------------------------------------------------------------------------#
library(ggplot2)
library(nlme)
library(multcomp)
library(dplyr)
#------------------------------------------------------------------------#
path.google <- "~/Google Drive/My Drive/MANDIFORE/MANDIFORE_CaseStudy_MortonArb/"

path.google <- "G:/.shortcut-targets-by-id/0B_Fbr697pd36c1dvYXJ0VjNPVms/MANDIFORE/MANDIFORE_CaseStudy_MortonArb/"

path.figures <- file.path(path.google, "Drought and heat analysis/Figures/SEA figures/")

runs.yr <- read.csv(file.path(path.google, "processed_data/All_runs_yearly.csv"))
runs.yr$Management <- car::recode(runs.yr$Management, "'None'='None'; 'Gap'='Group'; 'Shelter'='Shelter'; 'Under'='Under'")
runs.yr$Management <- factor(runs.yr$Management, levels = c("None", "Under", "Shelter", "Group"))
runs.yr$RCP.name <- car::recode(runs.yr$rcp, "'rcp45'='Low Emmissions'; 'rcp85'='High Emissions'")
runs.yr$RCP.name <- factor(runs.yr$RCP.name, levels=c("Low Emmissions", "High Emissions"))
runs.yr$loss.event.20 <- ifelse(runs.yr$agb.rel.diff<=-0.2, 1, 0)
summary(runs.yr)
head(runs.yr)

#Counting individual instances of a crash beginning
# CR Note: this was a hard-coded 5:nrow(), but I don't htink that was necessary nor functional since it only skipped the first entries for 1 GCM/RCP combo
for(i in 1:nrow(runs.yr)){ 
  GCM <- runs.yr[i, "GCM"]
  RCP <- runs.yr[i, "rcp"]
  MNG <- runs.yr[i, "Management"]
  YR <- runs.yr[i, "year"]
  if(YR > 2025){ # Need to skip 2025 because this code would reference 2025 off of 2024, which was a management year; this would mess with our analyses
    prev.20 <- runs.yr[runs.yr$GCM == GCM & runs.yr$rcp == RCP & runs.yr$Management == MNG & runs.yr$year == YR-1 , "loss.event.20"]
    runs.yr[i, "nonseq.loss.event.20"] <- ifelse((runs.yr[i, "loss.event.20"] == 1 & prev.20 ==F), 1, 0)
  }
}

runs.yr$crash <- ifelse(runs.yr$nonseq.loss.event.20==T, 1, 0)
summary(runs.yr)
runs.yr$density.tree.convert <- runs.yr$density.tree * 10000
#runs.yr <- runs.yr[runs.yr$year>=2025,]

#---------------------------------------------------#
# Here is why I start converting some of Christy's old script for our purposes
# This script pulls out duplicate situations into a different data.frame and adds them back at thend
# Duplicates refer to when one year can be a -1 year lag for one drought and a -5 for another
# group.crash.lag = time lag for GROUP of conditions with at least ONE RUN crashing (This was previously group.crash.lag)
# group.crash.lag.check --> Y/N indicating which set actually crashed (This was previously group.crash.lag.check)
# ind.crash.lag = time lag for individual management condition which crashed
#---------------------------------------------------#
# Extreme crash years
dups.df <- data.frame()
for(RCP in unique(runs.yr$rcp)){
  for(GCM in unique(runs.yr$GCM)){
    # Setting up some dummy columns
    runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM, "group.crash.lag"] <- NA
    runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM, "group.crash.lag.check"] <- "N"
    
    # For each RCP/GCM, checking for crash in each managment
    for(MNG in unique(runs.yr$Management)){
      runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG, "ind.crash.lag"] <- NA
      crash.event <- unique(runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$crash==1, "year"])
      crash.event <- sort(crash.event)
      for(YR in crash.event){
        
        temp.df <- data.frame()
        #Flagging when only one of the management scenarios has a crash
        if(!is.na(runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-5), "ind.crash.lag"])){
          temp.store <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-5),]
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-5), "ind.crash.lag"] <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-5), "year"] - YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-5), "crash.year"] <- YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-5), "group.crash.lag.check"] <- "Y"
          temp.df <- rbind(temp.df, temp.store)
        } else{
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-5), "ind.crash.lag"] <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-5), "year"] - YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-5), "crash.year"] <- YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-5), "group.crash.lag.check"] <- "Y"
        }
        
        if(!is.na(runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-4), "ind.crash.lag"])){
          temp.store <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-4),]
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-4), "ind.crash.lag"] <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-4), "year"] - YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-4), "crash.year"] <- YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-4), "group.crash.lag.check"] <- "Y"
          
          temp.df <- rbind(temp.df, temp.store)
        } else{
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-4), "ind.crash.lag"] <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-4), "year"] - YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-4), "crash.year"] <- YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-4), "group.crash.lag.check"] <- "Y"
        }
        
        if(!is.na(runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-3), "ind.crash.lag"])){
          temp.store <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-3),]
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-3), "ind.crash.lag"] <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-3), "year"] - YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-3), "crash.year"] <- YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-3), "group.crash.lag.check"] <- "Y"
          
          temp.df <- rbind(temp.df, temp.store)
        } else{
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-3), "ind.crash.lag"] <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-3), "year"] - YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-3), "crash.year"] <- YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-3), "group.crash.lag.check"] <- "Y"
        }
        
        if(!is.na(runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-2), "ind.crash.lag"])){
          temp.store <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-2),]
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-2), "ind.crash.lag"] <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-2), "year"] - YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-2), "crash.year"] <- YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-2), "group.crash.lag.check"] <- "Y"
          
          temp.df <- rbind(temp.df, temp.store)
        } else{
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-2), "ind.crash.lag"] <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-2), "year"] - YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-2), "crash.year"] <- YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-2), "group.crash.lag.check"] <- "Y"
        }
        
        if(!is.na(runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-1), "ind.crash.lag"])){
          temp.store <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-1),]
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-1), "ind.crash.lag"] <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-1), "year"] - YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-1), "crash.year"] <- YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR-1), "group.crash.lag.check"] <- "Y"
          
          temp.df <- rbind(temp.df, temp.store)
        } else{
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-1), "ind.crash.lag"] <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-1), "year"] - YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-1), "crash.year"] <- YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR-1), "group.crash.lag.check"] <- "Y"
        }
        
        if(!is.na(runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR), "ind.crash.lag"])){
          temp.store <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR),]
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR), "ind.crash.lag"] <- "loss"
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR), "crash.year"] <- YR
          temp.store[temp.store$rcp == RCP & temp.store$GCM == GCM & temp.store$Management == MNG & temp.store$year == (YR), "group.crash.lag.check"] <- "Y"
          temp.df <- rbind(temp.df, temp.store)
        } else{
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR), "ind.crash.lag"] <- "loss"
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR), "crash.year"] <- YR
          runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$Management == MNG & runs.yr$year == (YR), "group.crash.lag.check"] <- "Y"
        }
        
        dups.df <- rbind(dups.df, temp.df)
      }
    }
  }
}
summary(dups.df)
runs.yr <- runs.yr[runs.yr$year>=2025,]

#Adding the duplicate values back into the data frame
runs.yr <- rbind(runs.yr, dups.df)
runs.yr$ind.crash.lag <- factor(runs.yr$ind.crash.lag, levels=c("-5", "-4", "-3", "-2", "-1", "loss"))
summary(runs.yr)

runs.fill <- data.frame()
double <- 0
for(RCP in unique(runs.yr$rcp)){
  for(GCM in unique(runs.yr$GCM)){
    #runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM, "group.crash.lag"] <- NA
    crash.event <- unique(runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$crash==1, "year"])
    crash.event <- sort(crash.event)
    for(YR in crash.event){
      
      temp.df <- runs.yr[runs.yr$rcp == RCP & runs.yr$GCM == GCM & runs.yr$year >= (YR-5) & runs.yr$year <= (YR),]
      temp.df$ind.crash.lag <- ifelse((!is.na(temp.df$crash.year) & temp.df$crash.year != YR), NA, as.character(temp.df$ind.crash.lag))
      temp.df$group.crash.lag.check <- ifelse((!is.na(temp.df$crash.year) & temp.df$crash.year != YR), "N", temp.df$group.crash.lag.check)
      
      if(nrow(temp.df[temp.df$group.crash.lag.check == "Y",]) > 6){
        double <- double + 1
      }
      
      temp.df[temp.df$rcp == RCP & temp.df$GCM == GCM & temp.df$year >= (YR-5) & temp.df$year <= (YR), "crash.year"] <- YR
      temp.df[temp.df$rcp == RCP & temp.df$GCM == GCM & temp.df$year == (YR-5), "group.crash.lag"] <- -5
      temp.df[temp.df$rcp == RCP & temp.df$GCM == GCM & temp.df$year == (YR-4), "group.crash.lag"] <- -4
      temp.df[temp.df$rcp == RCP & temp.df$GCM == GCM & temp.df$year == (YR-3), "group.crash.lag"] <- -3
      temp.df[temp.df$rcp == RCP & temp.df$GCM == GCM & temp.df$year == (YR-2), "group.crash.lag"] <- -2
      temp.df[temp.df$rcp == RCP & temp.df$GCM == GCM & temp.df$year == (YR-1), "group.crash.lag"] <- -1
      temp.df[temp.df$rcp == RCP & temp.df$GCM == GCM & temp.df$year == YR, "group.crash.lag"] <- "loss"
      
      runs.fill <- rbind(runs.fill, temp.df)
    }
  }
}

runs.fill$ind.crash.lag <- factor(runs.fill$ind.crash.lag, levels=c("-5", "-4", "-3", "-2", "-1", "loss"))
runs.fill$group.crash.lag <- factor(runs.fill$group.crash.lag, levels=c("-5", "-4", "-3", "-2", "-1", "loss"))
summary(runs.fill)

time.weath.agg <- aggregate(cbind(diff.tair, rel.precip, rel.VPD)~ind.crash.lag, data = runs.fill[!is.na(runs.fill$ind.crash.lag),], FUN = mean, na.action = NULL)
time.weath.agg[, c("diff.tair.sd", "rel.precip.sd", "rel.VPD.sd")] <- aggregate(cbind(diff.tair, rel.precip, rel.VPD)~ind.crash.lag, data = runs.fill[!is.na(runs.fill$ind.crash.lag),], FUN = sd, na.action = NULL)[, c("diff.tair", "rel.precip", "rel.VPD")]

write.csv(time.weath.agg, file.path(path.google, "processed_data/Time_by_relweather.csv"), row.names = F)


#-----------------------------------------------------#
# Here is where we start running the analysis to make figures
#-----------------------------------------------------#
#-----------------------------------------------------#
# Looking at relative weather before a crash
# ind.crash.lag = time lag for individual management which crashed
# We include the crash year for this evaluation because we are working with temperature
#-----------------------------------------------------#
summary(runs.fill)

theme.clean <- theme(axis.text = element_text(size=rel(2), face="bold", color="black"),
                     axis.title = element_text(size=rel(1.7), face="bold"),
                     legend.title=element_text(size=rel(2)),
                     legend.text=element_text(size=rel(2)),
                     panel.background = element_rect(fill=NA, color="black"),
                     panel.grid=element_blank(),
                     panel.spacing.x = unit(1, "lines"),
                     plot.margin = unit(c(1, 1, 1, 1), "lines"))

raw.met.tair <- ggplot(data=runs.fill[!is.na(runs.fill$ind.crash.lag),], aes(x=ind.crash.lag, y=diff.tair, group=Management), position=dodge) +
  geom_errorbar(aes(color=Management), stat="summary", fun.y="sd", size=0.75, alpha=0.75) +
  geom_line(aes(color=Management), stat="summary", fun="mean", size=2) +
  geom_point(aes(color=Management), stat="summary", fun="mean", size=4) +
  scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  #theme(text=element_text(size=17))+
  theme.clean+
  ylab("Diff in Temperature (C)") +
  guides(color=F)

raw.met.precip <- ggplot(data=runs.fill[!is.na(runs.fill$ind.crash.lag),], aes(x=ind.crash.lag, y=rel.precip, group=Management), position=dodge) +
  geom_errorbar(aes(color=Management), stat="summary", fun.y="sd", size=0.75, alpha=0.75) +
  geom_line(aes(color=Management), stat="summary", fun="mean", size=2) +
  geom_point(aes(color=Management), stat="summary", fun="mean", size=4) +
  scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  #theme(text=element_text(size=17))+
  theme.clean+
  ylab("Relative Total Precip (m)")

raw.met.vpd <- ggplot(data=runs.fill[!is.na(runs.fill$ind.crash.lag),], aes(x=ind.crash.lag, y=rel.VPD, group=Management), position=dodge) +
  geom_errorbar(aes(color=Management), stat="summary", fun.y="sd", size=0.75, alpha=0.75) +
  geom_line(aes(color=Management), stat="summary", fun="mean", size=2) +
  geom_point(aes(color=Management), stat="summary", fun="mean", size=4) +
  scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  #theme(text=element_text(size=17))+
  theme.clean+
  ylab("Relative VPD (kPa)") + guides(color=F)


png(file.path(path.figures, "SEA_RelWeather_TimeMgmt_RawDat.png"), width=12, height=8, units="in", res=220)
cowplot::plot_grid(raw.met.tair, raw.met.precip, raw.met.vpd, ncol=2, labels = c("A", "B", "C"), rel_widths = c(0.8, 1.2))
dev.off()

relmet.var <- c("rel.precip", "diff.tair", "rel.VPD")
df.lag.relmetxind <- data.frame()
df.ano.relmetxind <- data.frame()
df.time.relmetxind <- data.frame()
df.mgmt.relmetxind <- data.frame()
df.plot.relmetxind <- data.frame()
for(COL in relmet.var){
  
  mod.lag <- lm(eval(substitute(j ~ relevel(as.factor(ind.crash.lag), "-5")*relevel(Management, "None"), list(j = as.name(COL)))), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
  # summary(mod.lag)$tTable
  
  df.ano <- data.frame(anova(mod.lag))
  # output <- summary(mod.lag)
  df.ano$comp <- rownames(df.ano)
  df.ano$VAR <- COL
  rownames(df.ano) <- NULL
  
  df.ano.relmetxind <- rbind(df.ano.relmetxind, df.ano)
  # 
  # # This permutation of the analysis will let us get the random effects adjusted means & SEs for each time/mgmt group to make a clean figure; do NOT interpret the p-values as they'll just indicate difference from 0
  # # NOTE: The random effects make HUGE differences on the directionality of means (e.g. w/o R.E. SHelter has *higher* has diff.tair, but when including R.E. it's *lower*)
  # # --> this means that it's important to be showing the mixed model estiamtes, BUT the SE associated with each point are similar and very large and make it really hard to tease out the stat sig that we describe
  # mod.plot <- nlme::lme(eval(substitute(j ~ relevel(ind.crash.lag, "-5")*relevel(Management, "None")- relevel(ind.crash.lag, "-5") - relevel(Management, "None") -1, list(j = as.name(COL)))), random=list(rcp = ~1, GCM =~1), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
  # 
  # 
  # df.plot <- data.frame(summary(mod.plot)$tTable)
  # df.plot$time <- c("-5", "-4", "-3", "-2", "-1", "loss")
  # df.plot$Management <- rep(c("None", "Under", "Shelter", "Group"), each=6)
  # df.plot$VAR <- COL
  # rownames(df.plot) <- NULL
  # 
  # df.plot.relmetxind <- rbind(df.plot.relmetxind, df.plot)
  # 
  
  # Doing the univariate analysis & saving output
  # time to crash
  mod.time <- nlme::lme(eval(substitute(j ~ relevel(as.factor(ind.crash.lag), "-5"), list(j = as.name(COL)))), random=list(Management=~1), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
  df.time <- data.frame(summary(mod.time)$tTable)
  
  # mod.time <- lm(eval(substitute(j ~ relevel(as.factor(ind.crash.lag), "-5"), list(j = as.name(COL)))), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
  # df.time <- data.frame(summary(mod.time)$coefficients)
  
  df.time$comp <- c("-5", "-4", "-3", "-2", "-1", "loss")
  df.time$VAR <- COL
  rownames(df.time) <- NULL
  
  df.time.relmetxind <- rbind(df.time.relmetxind, df.time)
  
  # management
  mod.mgmt <- nlme::lme(eval(substitute(j ~ relevel(Management, "None"), list(j = as.name(COL)))), random=list(ind.crash.lag=~1), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
  df.mgmt <- data.frame(summary(mod.mgmt)$tTable)
  
  # mod.mgmt <- lm(eval(substitute(j ~ relevel(Management, "None"), list(j = as.name(COL)))), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
  # df.mgmt <- data.frame(summary(mod.mgmt)$coefficients)
  
  df.mgmt$comp <- c("None", "Under", "Shelter", "Group")
  df.mgmt$VAR <- COL
  rownames(df.mgmt) <- NULL
  
  df.mgmt.relmetxind <- rbind(df.mgmt.relmetxind, df.mgmt)
  
  # eff.lag <- nlme::lme(eval(substitute(j ~ relevel(as.factor(ind.crash.lag), "-5")*relevel(Management, "None"), list(j = as.name(COL)))), random=list(rcp = ~1, GCM =~1), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
  # mod.lag.eff <- summary(mod.lag)$tTable
}

df.ano.relmetxind
# summary(df.lag.relmetxind)

# df.ano.relmetxind <- df.ano.relmetxind[,c(6,5,1,2,3,4)]
names(df.ano.relmetxind) <- car::recode(names(df.ano.relmetxind), c("'Pr..F.'='p.value'"))
df.ano.relmetxind <- droplevels(df.ano.relmetxind[!df.ano.relmetxind$comp=="Residuals",])
df.ano.relmetxind <- df.ano.relmetxind[,c("VAR", "comp", "Df", "F.value", "p.value")]
df.ano.relmetxind$comp <- gsub("relevel", "", df.ano.relmetxind$comp)
df.ano.relmetxind$comp <- gsub('"', "", df.ano.relmetxind$comp)
df.ano.relmetxind$comp <- gsub('[(]', "", df.ano.relmetxind$comp)
df.ano.relmetxind$comp <- gsub('[)]', "", df.ano.relmetxind$comp)
df.ano.relmetxind$comp <- gsub("as.factorind.crash.lag, -5", "Time", df.ano.relmetxind$comp)
df.ano.relmetxind$comp <- gsub("Management, None", "Harvest Scenario", df.ano.relmetxind$comp)
#IMPORTANT!!!!! Rounding for easy reading. This should not be how the values are reported in the end
df.ano.relmetxind$p.value <- round(df.ano.relmetxind$p.value, 3)
# View(df.ano.relmetxind)
df.ano.relmetxind

# names(df.time.relmetxind) <- car::recode(names(df.ano.relmetxind), c("'Pr...t..'='p.value'"))
df.time.relmetxind$p.value <- round(df.time.relmetxind$p.value, 3)
df.time.relmetxind

# names(df.mgmt.relmetxind) <- car::recode(names(df.mgmt.relmetxind), c("'Pr...t..'='p.value'"))
df.mgmt.relmetxind$p.value <- round(df.mgmt.relmetxind$p.value, 3)
df.mgmt.relmetxind

write.csv(df.ano.relmetxind, file.path(path.google, "Drought and heat analysis", "Mixed effects models results/SEA_ANOVA_RelMet_TimeMgmt.csv"), row.names = F)
write.csv(df.time.relmetxind, file.path(path.google, "Drought and heat analysis", "Mixed effects models results/SEA_ANOVA_RelMet_Time.csv"), row.names = F)
write.csv(df.mgmt.relmetxind, file.path(path.google, "Drought and heat analysis", "Mixed effects models results/SEA_ANOVA_RelMet_Mgmt.csv"), row.names = F)


# Extra analysis for tair harvest scenario based on results --> we need the multiple comparsions
tair.mgmt <- nlme::lme(diff.tair ~ Management, random=list(ind.crash.lag=~1), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
summary(tair.mgmt)$tTable

# Doing the actual pairwise comparisons
tair.mc <- multcomp::glht(tair.mgmt, linfct = mcp(Management = 'Tukey'))
summary(tair.mc)






# df.plot.relmetxind$time <- factor(df.plot.relmetxind$time, levels=c("-5", "-4", "-3", "-2", "-1", "loss"))
# ggplot(data=df.plot.relmetxind, aes(x=time, y=Value, group=Management, color=Management)) +
#   facet_wrap(~VAR, scales="free_y", ncol=2) +
#   geom_line(size=1.5) + 
#   geom_point(size=2) +
#   geom_errorbar(aes(x=time, ymin=Value-Std.Error, ymax=Value+Std.Error, group=Management, color=Management), size=1.5, alpha=0.5) +
#   scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
#   theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))






#-----------------------------------------------------------#
# Structural variables including those that didn't crash
#-----------------------------------------------------------#
#------#
# Just looking at the difference between those that did and didn't crash
#------#
raw.struc.agb <- ggplot(data=runs.fill[!is.na(runs.fill$group.crash.lag),], aes(x=group.crash.lag, y=agb, group=group.crash.lag.check), position=dodge) + 
  geom_rect(xmin=5.5, xmax=6.5, ymin=-Inf, ymax=Inf, fill="gray90", alpha=0.9, color=NA) + # I don't know why the alpha isn't working, but :shrug:
  coord_cartesian(ylim=c(5.5, 11.25)) +
  geom_errorbar(aes(color=group.crash.lag.check), stat="summary", fun.y="sd", size=1.25, alpha=0.5) +
  geom_line(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=1.5) +
  geom_point(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=2) +
  #scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  scale_color_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  theme(text=element_text(size=21))+
  ylab("AGB (kgC/m2)")+ guides(color=F)

raw.struc.density <- ggplot(data=runs.fill[!is.na(runs.fill$group.crash.lag),], aes(x=group.crash.lag, y=density.tree.convert, group=group.crash.lag.check), position=dodge) +
  geom_rect(xmin=5.5, xmax=6.5, ymin=-Inf, ymax=Inf, fill="gray90", alpha=0.9, color=NA) + # I don't know why the alpha isn't working, but :shrug:
  coord_cartesian(ylim=c(125, 300)) +
  geom_errorbar(aes(color=group.crash.lag.check), stat="summary", fun.y="sd", size=1.25, alpha=0.5) +
  geom_line(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=1.5) +
  geom_point(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=2) +
  #scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  scale_color_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  theme(text=element_text(size=21))+
  ylab("Density (trees/ha)")

raw.struc.meandbh <- ggplot(data=runs.fill[!is.na(runs.fill$group.crash.lag),], aes(x=group.crash.lag, y=tree.dbh.mean, group=group.crash.lag.check), position=dodge) +
  geom_rect(xmin=5.5, xmax=6.5, ymin=-Inf, ymax=Inf, fill="gray90", alpha=0.9, color=NA) + # I don't know why the alpha isn't working, but :shrug:
  coord_cartesian(ylim=c(31, 35)) +
  geom_errorbar(aes(color=group.crash.lag.check), stat="summary", fun.y="sd", size=1.25, alpha=0.5) +
  geom_line(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=1.5) +
  geom_point(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=2) +
  #scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  scale_color_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  theme(text=element_text(size=21))+
  ylab("Mean DBH (cm)")+ guides(color=F)

raw.struc.sddbh <- ggplot(data=runs.fill[!is.na(runs.fill$group.crash.lag),], aes(x=group.crash.lag, y=tree.dbh.sd, group=group.crash.lag.check), position=dodge) +
  geom_rect(xmin=5.5, xmax=6.5, ymin=-Inf, ymax=Inf, fill="gray90", alpha=0.9, color=NA) + # I don't know why the alpha isn't working, but :shrug:
  coord_cartesian(ylim=c(16.75, 19.25)) +
  geom_errorbar(aes(color=group.crash.lag.check), stat="summary", fun.y="sd", size=1.25, alpha=0.5) +
  geom_line(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=1.5) +
  geom_point(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=2) +
  #scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  scale_color_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  theme(text=element_text(size=21))+
  ylab("SD of DBH (cm)")


png(file.path(path.figures, "SEA_Structure_TimeCrashYN_RawDat.png"), width=12, height=8, units="in", res=220)
cowplot::plot_grid(raw.struc.agb, raw.struc.density, raw.struc.meandbh, raw.struc.sddbh, ncol=2, rel_widths = c(1,1.5,1,1.5))
dev.off()









#---------------------------------#
# Adding Management as a facet wrapped value
#---------------------------------#

raw.struc.agb2 <- ggplot(data=runs.fill[!is.na(runs.fill$group.crash.lag),], aes(x=group.crash.lag, y=agb, group=group.crash.lag.check), position=dodge) +
  geom_rect(xmin=5.5, xmax=6.5, ymin=-Inf, ymax=Inf, fill="gray90", alpha=0.9, color=NA) + # I don't know why the alpha isn't working, but :shrug:
  coord_cartesian(ylim=c(3.8, 12.5)) +
  facet_wrap(~Management)+
  geom_errorbar(aes(color=group.crash.lag.check), stat="summary", fun.y="sd", size=1.25, alpha=0.5) +
  geom_line(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=1.5) +
  geom_point(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=2) +
  #scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  scale_color_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  ylab("AGB (kgC/m2)")

png(file.path(path.figures, "SEA_Structure-AGB_TimeCrashYN_RawDat.png"), width=8, height=8, units="in", res=220)
raw.struc.agb2
dev.off()


raw.struc.density2 <- ggplot(data=runs.fill[!is.na(runs.fill$group.crash.lag),], aes(x=group.crash.lag, y=density.tree.convert, group=group.crash.lag.check), position=dodge) +
  geom_rect(xmin=5.5, xmax=6.5, ymin=-Inf, ymax=Inf, fill="gray90", alpha=0.9, color=NA) + # I don't know why the alpha isn't working, but :shrug:
  coord_cartesian(ylim=c(100, 390)) +
  facet_wrap(~Management)+
  geom_errorbar(aes(color=group.crash.lag.check), stat="summary", fun.y="sd", size=1.25, alpha=0.5) +
  geom_line(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=1.5) +
  geom_point(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=2) +
  #scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  scale_color_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  ylab("Density (trees/ha)")

png(file.path(path.figures, "SEA_Structure-Density_TimeCrashYN_RawDat.png"), width=8, height=8, units="in", res=220)
raw.struc.density2
dev.off()


raw.struc.meandbh2 <- ggplot(data=runs.fill[!is.na(runs.fill$group.crash.lag),], aes(x=group.crash.lag, y=tree.dbh.mean, group=group.crash.lag.check), position=dodge) +
  geom_rect(xmin=5.5, xmax=6.5, ymin=-Inf, ymax=Inf, fill="gray90", alpha=0.9, color=NA) + # I don't know why the alpha isn't working, but :shrug:
  coord_cartesian(ylim=c(24, 40)) +
  facet_wrap(~Management)+
  geom_errorbar(aes(color=group.crash.lag.check), stat="summary", fun.y="sd", size=1.25, alpha=0.5) +
  geom_line(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=1.5) +
  geom_point(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=2) +
  #scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  scale_color_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  ylab("Mean DBH (cm)")

png(file.path(path.figures, "SEA_Structure-DBHmean_TimeCrashYN_RawDat.png"), width=8, height=8, units="in", res=220)
raw.struc.meandbh2
dev.off()

raw.struc.sddbh2 <- ggplot(data=runs.fill[!is.na(runs.fill$group.crash.lag),], aes(x=group.crash.lag, y=tree.dbh.sd, group=group.crash.lag.check), position=dodge) +
  geom_rect(xmin=5.5, xmax=6.5, ymin=-Inf, ymax=Inf, fill="gray90", alpha=0.9, color=NA) + # I don't know why the alpha isn't working, but :shrug:
  coord_cartesian(ylim=c(10, 22)) +
  facet_wrap(~Management)+
  geom_errorbar(aes(color=group.crash.lag.check), stat="summary", fun.y="sd", size=1.25, alpha=0.5) +
  geom_line(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=1.5) +
  geom_point(aes(color=group.crash.lag.check), stat="summary", fun="mean", size=2) +
  #scale_color_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  scale_color_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  ylab("SD of DBH (cm)")

png(file.path(path.figures, "SEA_Structure-DBHsd_TimeCrashYN_RawDat.png"), width=8, height=8, units="in", res=220)
raw.struc.sddbh2
dev.off()

plot(runs.fill$tree.dbh.sd[!is.na(runs.fill$group.crash.lag)] ~ runs.fill$tree.dbh.mean[!is.na(runs.fill$group.crash.lag)])


struct.var <- c("agb", "density.tree.convert", "tree.dbh.mean", "tree.dbh.sd")
df.lag.structxind <- data.frame()
df.ano.structxind <- data.frame()
df.time.structxind <- data.frame()
df.mgmt.structxind <- data.frame()
df.plot.structxind <- data.frame()

runs.fill$group.crash.lag.check <- as.factor(runs.fill$group.crash.lag.check)
for(COL in struct.var){
  # #   mod.lag <- nlme::lme(eval(substitute(j ~ relevel(as.factor(group.crash.lag), "-1")*relevel(as.factor(group.crash.lag.check), "N")*relevel(Management, "None"), list(j = as.name(COL)))), random=list(rcp = ~1, GCM =~1), data = runs.fill[!is.na(runs.fill$group.crash.lag) & runs.fill$group.crash.lag!="loss",], na.action = na.omit)
  # #   anova(mod.lag)
  
  mod.lag <- lm(eval(substitute(j ~ relevel(as.factor(group.crash.lag), "-5")*relevel(Management, "None")*relevel(as.factor(group.crash.lag.check), "N"), list(j = as.name(COL)))), data = runs.fill[runs.fill$group.crash.lag!="loss",], na.action = na.omit)
  # summary(mod.lag)$tTable
  
  df.ano <- data.frame(anova(mod.lag))
  # output <- summary(mod.lag)
  df.ano$comp <- rownames(df.ano)
  df.ano$VAR <- COL
  rownames(df.ano) <- NULL
  
  df.ano.structxind <- rbind(df.ano.structxind, df.ano)
  # 
  
  # Doing the univariate analysis & saving output
  # time to crash
  mod.time <- nlme::lme(eval(substitute(j ~ relevel(as.factor(group.crash.lag), "-5")*relevel(as.factor(group.crash.lag.check), "N"), list(j = as.name(COL)))), random=list(Management=~1), data = runs.fill[runs.fill$group.crash.lag!="loss",], na.action = na.omit)
  df.time <- data.frame(summary(mod.time)$tTable)
  
  # mod.time <- lm(eval(substitute(j ~ relevel(as.factor(ind.crash.lag), "-5"), list(j = as.name(COL)))), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
  # df.time <- data.frame(summary(mod.time)$coefficients)
  
  df.time$time <- c("-5", "-4", "-3", "-2", "-1")
  df.time$crash <- rep(c("N", "Y"), each=5)
  df.time$VAR <- COL
  rownames(df.time) <- NULL
  
  df.time.structxind <- rbind(df.time.structxind, df.time)
  
  # management
  mod.mgmt <- nlme::lme(eval(substitute(j ~ relevel(Management, "None")*relevel(as.factor(group.crash.lag.check), "N")-1, list(j = as.name(COL)))), random=list(group.crash.lag=~1), data = runs.fill[runs.fill$group.crash.lag!="loss",], na.action = na.omit)
  df.mgmt <- data.frame(summary(mod.mgmt)$tTable)
  
  # mod.mgmt <- lm(eval(substitute(j ~ relevel(Management, "None"), list(j = as.name(COL)))), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
  # df.mgmt <- data.frame(summary(mod.mgmt)$coefficients)
  
  df.mgmt$mgmt <- c("None", "Under", "Shelter", "Group") # Caution: this is hard-coded!
  df.mgmt$crash <- rep(c("N", "Y"), each=4)
  df.mgmt$VAR <- COL
  rownames(df.mgmt) <- NULL
  
  df.mgmt.structxind <- rbind(df.mgmt.structxind, df.mgmt)
  
  # eff.lag <- nlme::lme(eval(substitute(j ~ relevel(as.factor(ind.crash.lag), "-5")*relevel(Management, "None"), list(j = as.name(COL)))), random=list(rcp = ~1, GCM =~1), data = runs.fill[!is.na(runs.fill$ind.crash.lag),], na.action = na.omit)
  # mod.lag.eff <- summary(mod.lag)$tTable
}

df.ano.structxind
# summary(df.lag.structxind)

# df.ano.structxind <- df.ano.structxind[,c(6,5,1,2,3,4)]
names(df.ano.structxind) <- car::recode(names(df.ano.structxind), c("'Pr..F.'='p.value'"))
df.ano.structxind <- droplevels(df.ano.structxind[!df.ano.structxind$comp=="Residuals",])
df.ano.structxind <- df.ano.structxind[,c("VAR", "comp", "Df", "F.value", "p.value")]
df.ano.structxind$comp <- gsub("relevel", "", df.ano.structxind$comp)
df.ano.structxind$comp <- gsub('"', "", df.ano.structxind$comp)
df.ano.structxind$comp <- gsub('[(]', "", df.ano.structxind$comp)
df.ano.structxind$comp <- gsub('[)]', "", df.ano.structxind$comp)
df.ano.structxind$comp <- gsub("as.factorgroup.crash.lag, -5", "Time", df.ano.structxind$comp)
df.ano.structxind$comp <- gsub("as.factorgroup.crash.lag.check, N", "Crash", df.ano.structxind$comp)
df.ano.structxind$comp <- gsub("Management, None", "Harvest Scenario", df.ano.structxind$comp)
#IMPORTANT!!!!! Rounding for easy reading. This should not be how the values are reported in the end
df.ano.structxind$p.value <- round(df.ano.structxind$p.value, 3)
# View(df.ano.structxind)
df.ano.structxind

# names(df.time.structxind) <- car::recode(names(df.ano.structxind), c("'Pr...t..'='p.value'"))
df.time.structxind$p.value <- round(df.time.structxind$p.value, 3)
df.time.structxind

# names(df.mgmt.structxind) <- car::recode(names(df.mgmt.structxind), c("'Pr...t..'='p.value'"))
df.mgmt.structxind$p.value <- round(df.mgmt.structxind$p.value, 3)
df.mgmt.structxind

write.csv(df.ano.structxind, file.path(path.google, "Drought and heat analysis", "Mixed effects models results/SEA_ANOVA_Struct_TimeMgmt.csv"), row.names = F)
write.csv(df.time.structxind, file.path(path.google, "Drought and heat analysis", "Mixed effects models results/SEA_ANOVA_Struct_Time.csv"), row.names = F)
write.csv(df.mgmt.structxind, file.path(path.google, "Drought and heat analysis", "Mixed effects models results/SEA_ANOVA_Struct_Mgmt.csv"), row.names = F)


theme.clean <- theme(axis.text = element_text(size=rel(1.7), face="bold", color="black"),
                     axis.title = element_text(size=rel(1.7), face="bold"),
                     legend.title=element_text(size=rel(2)),
                     legend.text=element_text(size=rel(2)),
                     panel.background = element_rect(fill=NA, color="black"),
                     panel.grid=element_blank(),
                     panel.spacing.x = unit(1, "lines"),
                     plot.margin = unit(c(1, 1, 1, 1), "lines"))

raw.struc.agb3 <- ggplot(data=runs.fill[runs.fill$group.crash.lag!="loss",], aes(x=Management, fill=group.crash.lag.check, y=agb), position=dodge) + 
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  #theme(text=element_text(size=21))+
  ylab("AGB (kgC/m2)")+
  theme.clean+
  guides(fill=F)

raw.struc.density3 <- ggplot(data=runs.fill[runs.fill$group.crash.lag!="loss",], aes(x=Management, fill=group.crash.lag.check, y=density.tree.convert), position=dodge) + 
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+
  #theme(text=element_text(size=21))+
  theme.clean+
  ylab("Density (trees/ha)")

raw.struc.meandbh3 <- ggplot(data=runs.fill[runs.fill$group.crash.lag!="loss",], aes(x=Management, fill=group.crash.lag.check, y=tree.dbh.mean), position=dodge) + 
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+  
  #theme(text=element_text(size=21))+
  theme.clean+
  ylab("Mean DBH (cm)")+
  guides(fill=F)

raw.struc.sddbh3 <- ggplot(data=runs.fill[runs.fill$group.crash.lag!="loss",], aes(x=Management, fill=group.crash.lag.check, y=tree.dbh.sd), position=dodge) + 
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(name = "Loss Event\nOccurence", values=c("Y"="red3", "N"="black")) +
  theme_bw() + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))+ 
  theme.clean+
  #theme(text=element_text(size=21))+
  ylab("SD of DBH (cm)")


png(file.path(path.figures, "SEA_Structure_CrashYN_RawDat-Boxplot.png"), width=12, height=8, units="in", res=220)
  cowplot::plot_grid(raw.struc.agb3, raw.struc.density3, raw.struc.meandbh3, raw.struc.sddbh3, ncol=2, labels = c("A", "B", "C", "D"), rel_widths = c(0.8,1.2,0.8,1.2))
dev.off()

## Summarizing results
df.mgmt.structxind[order(df.mgmt.structxind$mgmt),]

mod.loss <- nlme::lme(density.tree.convert ~ relevel(as.factor(group.crash.lag.check), "N"), random=list(Management = ~1, group.crash.lag=~1), data = runs.fill[runs.fill$group.crash.lag!="loss",], na.action = na.omit)
summary(mod.loss)$tTable

mod.mgmt <- nlme::lme(density.tree.convert ~ Management, random=list(group.crash.lag.check = ~1, group.crash.lag=~1), data = runs.fill[runs.fill$group.crash.lag!="loss",], na.action = na.omit)
summary(mod.mgmt)$tTable

mgmt.mc <- multcomp::glht(mod.mgmt, linfct = mcp(Management = 'Tukey'))
summary(mgmt.mc)