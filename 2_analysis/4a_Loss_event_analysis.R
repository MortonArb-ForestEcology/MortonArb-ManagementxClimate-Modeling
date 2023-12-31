#----------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick, Christy Rollinson
# Project: Mandifore -- Morton Arb Management Case Study
# Purpose: This script runs AIC model selection, does linear regression, and creates figures
# Inputs: Yearly ED2 output csv from script 2a_Yearly_ED2_sum.R
# Outputs: Figures and Tables
# Notes: 
#----------------------------------------------------------------------------------------------------------------------#
library(lubridate)
library(dplyr)
library(lubridate)
library(nlme)
library(AICcmodavg)
library(ggplot2)
library(ggpubr)
library(multcomp)

path.google <- "G:/.shortcut-targets-by-id/0B_Fbr697pd36c1dvYXJ0VjNPVms/MANDIFORE/MANDIFORE_CaseStudy_MortonArb/"
path.figures <- file.path(path.google, "Drought and heat analysis/Figures/")

runs.comb <- read.csv(paste0(path.google, "processed_data/All_runs_yearly.csv"))

runs.comb$Management <- car::recode(runs.comb$Management, "'None'='None'; 'Gap'='Group'; 'Shelter'='Shelter'; 'Under'='Under'")
runs.comb$Management <- factor(runs.comb$Management, levels = c("None", "Group", "Shelter", "Under"))

runs.comb$Driver.set <- paste0(runs.comb$GCM,"." ,runs.comb$rcp)

runs.comb$RCP.name <- car::recode(runs.comb$rcp, "'rcp45'='RCP4.5'; 'rcp85'='RCP8.5'")
runs.comb$RCP.name <- factor(runs.comb$RCP.name, levels=c("RCP4.5", "RCP8.5"))

runs.late <- runs.comb[runs.comb$year >= 2025,]

runs.late <- runs.late[!is.na(runs.late$agb.rel.diff.future),]

#-------------------------------------------------------------#
#Counting the number of major agb loss events
#-------------------------------------------------------------#
#agb.extreme <- as.numeric(quantile(runs.late$agb.rel.diff.future, probs = c(.025)))

runs.comb$loss.event.20 <- ifelse(runs.comb$agb.rel.diff <= -.20, 1, 0)

#Counting individual instances of a crash beginning
for(i in 5:nrow(runs.comb)){
  DRIVE <- runs.comb[i, "Driver.set"]
  MNG <- runs.comb[i, "Management"]
  YR <- runs.comb[i, "year"]
  if(YR != 2007){
    prev.20 <- runs.comb[runs.comb$Driver.set == DRIVE & runs.comb$Management == MNG & runs.comb$year == YR-1 , "loss.event.20"]
    runs.comb[i, "nonseq.loss.event.20"] <- ifelse((runs.comb[i, "loss.event.20"] == 1 & prev.20 ==F), 1, 0)
  }
}

#Modified so that red dots are now looking at unique occurences of a crash instead of the total number of years crashing
#--------------------------------------#
# Figure 4
#--------------------------------------#
levels(runs.comb$Management) <- c("None", "Under", "Shelter", "Group")
png(paste0(path.figures, "AGB_static_crashes.png"), width=16, height=8, units="in", res=220)
ggplot(data=runs.comb)+
  facet_grid(Management ~ RCP.name) +
  geom_rect(xmin=2020, xmax=2024, ymin=-Inf, ymax=Inf, fill="gray73", alpha=0.9) +
  geom_line(aes(x=year, y=agb, group=GCM)) +
  geom_vline(aes(xintercept=2050), linetype = "dashed")+
  geom_vline(aes(xintercept=2099), linetype = "dashed")+
  #geom_point(data=runs.comb[runs.comb$loss.event.20==T & runs.comb$year>2024 & !is.na(runs.comb$agb.rel.diff),], aes(x=year, y=agb, group=GCM), color="blue2", size=3) +
  geom_point(data=runs.comb[runs.comb$nonseq.loss.event.20==T & runs.comb$year>2024 & !is.na(runs.comb$agb.rel.diff),], aes(x=year, y=agb, group=GCM), color="red2", size=3) +
  labs(x="Year", y="Aboveground Biomass (kgC/m2)") +
  #geom_text(x=2025, y=25, label="Harvest Period: 2020-2024", color="orange3", hjust=0) +
  theme(axis.text = element_text(size=rel(2), color="black"),
        axis.title = element_text(size=rel(2), face="bold"),
        panel.background = element_rect(fill=NA, color="black"),
        panel.grid=element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.text = element_text(size=rel(2), face="bold"),
        strip.background = element_rect(fill=NA),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        plot.title = element_text(size=rel(2), face="bold", hjust=0.5))
dev.off()

# Getting summaries of # crashes by scenario by mid & end of century
runs.comb$crash <- ifelse(runs.comb$nonseq.loss.event.20==T, 1, 0)

crash.summary <- aggregate(crash~ GCM + rcp + RCP.name + Management, data=runs.comb[runs.comb$year>=2025 ,], FUN=sum)
names(crash.summary)[names(crash.summary)=="crash"] <- "crash.end"
crash.summary$crash.mid <- aggregate(crash~ GCM + rcp + RCP.name + Management, data=runs.comb[runs.comb$year>=2025 & runs.comb$year<=2050,], FUN=sum)$crash
summary(crash.summary)

library(MASS)
library(lme4)
crash.end <- glm.nb(crash.end~Management*rcp, data = crash.summary)
summary(crash.end)
anova(crash.end)

t.test(crash.summary$crash.mid[crash.summary$Management=="None"], crash.summary$crash.mid[crash.summary$Management=="Shelter"], paired=T)
t.test(crash.summary$crash.mid[crash.summary$Management=="Group"], crash.summary$crash.mid[crash.summary$Management=="Shelter"], paired=T)
t.test(crash.summary$crash.mid[crash.summary$Management=="Under"], crash.summary$crash.mid[crash.summary$Management=="Shelter"], paired=T)

# t.test(crash.summary$crash.mid[crash.summary$Management=="None"], crash.summary$crash.mid[crash.summary$Management=="Group"], paired=T)

# TukeyHSD(crash.mid)
# chisq.test(x=crash.summary$crash.mid[crash.summary$Management=="None"], y=crash.summary$crash.mid[crash.summary$Management=="Shelter"])

# Doesn't converge
# crash.end.mm <- glmer.nb(crash.end~Management*rcp + (1|GCM), data = crash.summary)
# summary(crash.end.mm)
crash.mid <- glm.nb(crash.mid~Management*rcp, data = crash.summary)
summary(crash.mid)
anova(crash.mid)

crash.stack <- stack(crash.summary[,c("crash.mid", "crash.end")])
crash.stack[,c("Management", "rcp", "GCM", "RCP.name")] <- crash.summary[,c("Management", "rcp", "GCM", "RCP.name")]
crash.stack$TimePeriod <- car::recode(crash.stack$ind, "'crash.mid'='mid-century'; 'crash.end'='end-of-century'")
crash.stack$TimePeriod <- factor(crash.stack$TimePeriod, levels=c("mid-century", "end-of-century"))

crash.prop <- crash.stack
crash.prop$values <- ifelse(crash.prop$values>=1, 1, 0)
stat.df <- aggregate(values~TimePeriod+Management+rcp, data = crash.prop, FUN = sum)
stat.shape <- reshape(stat.df,idvar =c("Management", "TimePeriod"), timevar = "rcp",direction = "wide")
colnames(stat.shape) <- c("TimePeriod", "Management", "rcp45", "rcp85")
stat.shape <- reshape(stat.shape,idvar ="Management", timevar = "TimePeriod",direction = "wide")
for(i in 2:ncol(stat.shape)){
  stat.shape[,i] <- paste0(stat.shape[,i], "/14")
}

stat.df <- aggregate(values~TimePeriod+Management+rcp, data = crash.stack, FUN = sum)
stat.shape <- reshape(stat.df,idvar =c("Management", "TimePeriod"), timevar = "rcp",direction = "wide")
colnames(stat.shape) <- c("TimePeriod", "Management", "rcp45", "rcp85")
stat.shape <- reshape(stat.shape,idvar ="Management", timevar = "TimePeriod",direction = "wide")

# ggplot(data=crash.stack) +
# facet_grid(.~rcp) +
# geom_boxplot(aes(x=ind, y=values, fill=Management))
theme.clean <-   theme(axis.text = element_text(size=rel(2), color="black"),
                       axis.title = element_text(size=rel(2), face="bold"),
                       legend.title=element_text(size=rel(1.5)),
                       legend.text=element_text(size=rel(1.5)),
                       legend.position = "top",
                       panel.background = element_rect(fill=NA, color="black"),
                       panel.grid=element_blank(),
                       panel.spacing.x = unit(1, "lines"),
                       strip.text = element_text(size=rel(2), face="bold"),
                       strip.background = element_rect(fill=NA),
                       plot.margin = unit(c(1, 1, 1, 1), "lines"),
                       plot.title = element_text(size=rel(2), face="bold", hjust=0.5))
  
#--------------------------------------#
# Figure 5
#--------------------------------------#
levels(crash.stack$Management) <- c("None", "Under", "Shelter", "Group")
png(paste0(path.figures, "Crashes_Summary.png"), width=10, height=8, units="in", res=220)
ggplot(data=crash.stack) +
  facet_grid(TimePeriod~.) +
  geom_bar(aes(x=rcp, y=values, fill=Management), position="dodge", stat="summary", fun.y="mean") +
  geom_errorbar(aes(x=rcp, y=values, fill=Management), position="dodge", stat="summary", fun.y="sd") +
  scale_y_continuous(name="Mean Number of Loss Events", expand=c(0,0)) +
  scale_x_discrete(labels=c("RCP4.5", "RCP8.5")) +
  coord_cartesian(ylim=c(0,2.1)) +
  scale_fill_manual(values=c("None"="#1f78b4", "Under"="#a6cee3", "Shelter"="#33a02c", "Group"="#b2df8a")) +
  theme.clean + theme(axis.title.x=element_blank(), panel.spacing.y = unit(2, "lines"))
dev.off()
