# PACKAGES
library(lme4);library(tidyverse); library(dplyr);library(writexl);library(ggpubr);library(lmerTest)
library(tidyr);library(car);library(RColorBrewer);library(ggplot2) #library(plyr);

rm(list=ls())
dodge.posn <- position_dodge(.9)

# PAPERS
theme_angele_ss <- theme(panel.background = element_blank(),
                         panel.border =element_rect(colour="black", fill=NA),
                         plot.background = element_blank(),
                         panel.grid = element_blank(),
                         axis.line = element_line(colour ="black"),
                         axis.text.x = element_text (size = 15,colour= "black", family="sans"),
                         axis.text.y = element_text (size = 15,colour= "black", family="sans"),
                         axis.ticks.y = element_line(colour="black"),
                         axis.ticks.x = element_line(colour="black"),
                         axis.title.x = element_text(size = 15, vjust = -0.5, family="sans"),
                         axis.title.y = element_text(size = 15, vjust = 2, family="sans"),
                         legend.text=  element_text(size = 15, family="sans", margin = margin(t = 10)),
                         legend.title = element_text(size = 15, vjust = 2, family="sans"),
                         legend.key = element_blank(),
                         legend.position = "right",
                         legend.spacing.x = unit(0.2, 'cm'),
                         title = element_text(size = 20, family="sans"),
                         strip.text = element_text(size = 15),
)
_______________________________________________________________________________________________________________________________________________________________________________________________
# PREPARE DATA

Data$SexRec <- str_sub(Data$SexCombinaison, start = 2)
Data <- Data %>% unite("SexSetting", c(Sex,Setting), sep = "/", remove = FALSE)

# Take out the Keepers-directed signals, None and NA in Recipient
Data <- Data[Data$Recipient != "Keepers", ]; Data <- Data[Data$Recipient != "NA", ]; Data <- Data[Data$Recipient != "None", ]
Data <- subset(Data, !is.na(Recipient))
# Take out lines "no focal subject"
Data <- Data[Data$Subject != "No focal subject", ] 

# Take out "Special" in Behaviors
Data <- Data[Data$Behavior != "Special", ] ; Data <- Data[is.na(Data$Behavior2) | Data$Behavior2 != "Special", ]
Data$Behavior[Data$Behavior == "Pout"] <- "Pant hoot face" ; Data$Behavior2[Data$Behavior2 == "Pout"] <- "Pant hoot face" ; Data$Behavior3[Data$Behavior3 == "Pout"] <- "Pant hoot face" ; Data$Behavior4[Data$Behavior4 == "Pout"] <- "Pant hoot face" ; Data$Behavior5[Data$Behavior5 == "Pout"] <- "Pant hoot face"
Data$Behavior[Data$Behavior == "Present grooming"] <- "Present" ; Data$Behavior2[Data$Behavior2 == "Present grooming"] <- "Present" ; Data$Behavior3[Data$Behavior3 == "Present grooming"] <- "Present" ; Data$Behavior4[Data$Behavior4 == "Present grooming"] <- "Present" ; Data$Behavior5[Data$Behavior5 == "Present grooming"] <- "Present"
Data$Behavior[Data$Behavior == "Present sexual"] <- "Present" ; Data$Behavior2[Data$Behavior2 == "Present sexual"] <- "Present" ; Data$Behavior3[Data$Behavior3 == "Present sexual"] <- "Present" ; Data$Behavior4[Data$Behavior4 == "Present sexual"] <- "Present" ; Data$Behavior5[Data$Behavior5 == "Present sexual"] <- "Present" 
Data$Behavior[Data$Behavior == "Present climb on me"] <- "Present" ; Data$Behavior2[Data$Behavior2 == "Present climb on me"] <- "Present" ; Data$Behavior3[Data$Behavior3 == "Present climb on me"] <- "Present" ; Data$Behavior4[Data$Behavior4 == "Present climb on me"] <- "Present" ; Data$Behavior5[Data$Behavior5 == "Present climb on me"] <- "Present"


# Only keep behaviors that are used at least twice by each ind
dat.ind.bhv <- data.frame(table(Data$Subject,Data$Behavior)) ; colnames(dat.ind.bhv) <- c("Subject","Behavior","Nsignals")
dat.ind.bhv$Behavior <- as.character(dat.ind.bhv$Behavior)
dat.ind.bhv <- subset(dat.ind.bhv, Nsignals>=2) # Keeps only the values in "column" that appear more than 2 times (N >= 2)
Data <- merge(dat.ind.bhv, Data, by=c("Subject", "Behavior"))

# Omit all levels of Subject that contributed fewer than 30 cases --- 117 --> 93 individuals
nb.signals.per.ind <- data.frame(Data %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(nb.signals.per.ind) <- c("Subject","TotComActs") 
nb.signals.per.ind$Subject <- as.character(nb.signals.per.ind$Subject)
nb.signals.per.ind <- subset(nb.signals.per.ind, nb.signals.per.ind$TotComActs>=30) # Keeps only the values in "column" that appear more than 30 times (N >= 30)
Data <- subset(Data, Subject %in% nb.signals.per.ind$Subject)

Data <- Data[!is.na(Data$Group), ]


# Number of interaction partners
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)
# Get list of adults (from Subject column)
adults <- unique(Data$Subject)
# Build undirected dyads (order doesn't matter)
dyads <- Data %>%
  filter(!is.na(Subject) & !is.na(Recipient)) %>%
  mutate(Partner1 = pmin(as.character(Subject), as.character(Recipient)),
         Partner2 = pmax(as.character(Subject), as.character(Recipient))) %>%
  distinct(Partner1, Partner2)
# Reshape to a long format with Individual and their Partner
dyads_long <- dyads %>%
  mutate(ID1 = Partner1, ID2 = Partner2) %>%
  select(ID1, ID2) %>%
  dplyr::rename(Individual = ID1, Partner = ID2) %>%
  bind_rows(
    dyads %>%
      mutate(ID1 = Partner2, ID2 = Partner1) %>%
      select(ID1, ID2) %>%
      dplyr::rename(Individual = ID1, Partner = ID2)
  )
# Count unique partners for adults only
NbInteractionPartners <- dyads_long %>%
  filter(Individual %in% adults) %>%
  group_by(Individual) %>%
  summarise(NbInteractionPartners = n_distinct(Partner), .groups = "drop")
# Join back into your main Data
Data <- Data %>%
  left_join(NbInteractionPartners, by = c("Subject" = "Individual"))

# Sampling effort
Datareduit <- Data[Data$Duplicated != "T", ]  
SamplingEffort <- data.frame(Datareduit %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(SamplingEffort) <- c("Subject","TotComActs") 
Data <- merge(SamplingEffort, Data, by=c("Subject"))

# Individual ASO size 
Data <- Data[Data$ASO != "NA", ]; Data <- Data[!is.na(Data$ASO), ]
IndASOSize <- aggregate(data=Data, ASO ~ Subject, function(x) length(unique(x))) ; colnames(IndASOSize) <- c("Subject","ASOSize") 
Data <- merge(IndASOSize, Data, by=c("Subject"))

IndASOSizeGroup <- aggregate(data=Data, ASO ~ Group, function(x) length(unique(x))) ; colnames(IndASOSizeGroup) <- c("Group","ASOSize") 

length(unique(Data[["Behavior"]])) #69
Data %>%
  filter(!is.na(`Behavioural category`), !is.na(Behavior)) %>%
  group_by(`Behavioural category`) %>%
  summarise(
    Nb = n_distinct(Behavior),
    .groups = "drop"
  )
#  Behavioural category      Nb Type 
#1 Bodily gesture            16
#2 Facial expression          6
#3 Limb gesture              34  
#4 Vocalization              13


________________________________________________________________________________________________________________________________________________________________________________________________
1. Cumulative repertoire : Curve to check wether the entire group repertoires have been grasped

library(stringr)

behaviorsID <- aggregate(data=Data, Behavior ~ ObservationID, function(x) length(unique(x))) ; colnames(behaviorsID) <- c("ID","N")
RepertoireID <- distinct(Data,Behavior,ObservationID)

MediaFile <- Data[,c("Media file","Group","Behavior")]

# Leipzig A
Leipzig <- Data[,c("ObservationID","Group","Behavior")]
LeipzigA <- Leipzig[(Leipzig$Group=="A"),] 
LeipzigA$Date <- str_sub(LeipzigA$ObservationID, start = 1, end = 4)
AbehaviorsID <- aggregate(data=LeipzigA, Behavior ~ Date, function(x) length(unique(x))) ; colnames(AbehaviorsID) <- c("Date","N")
LeipzigA <- merge(LeipzigA, AbehaviorsID, by=c("Date"))
ARepertoireID <- distinct(LeipzigA,Behavior,Date)

require(data.table)
LeipzigA2 <- as.data.table(unique(ARepertoireID))
setkey(LeipzigA2, "Date")
LeipzigA2[, Behavior := as.numeric(factor(Behavior, levels = unique(Behavior)))]
setkey(LeipzigA2, "Date", "Behavior")
LeipzigA2.out <- LeipzigA2[J(unique(Date)), mult="last"]
LeipzigA2.out[, Behavior := cummax(Behavior)]

ggplot(LeipzigA2.out, aes(Date, Behavior)) + 
  geom_line(aes(group = 1)) +
  geom_point(size=3) + 
  theme_angele_ss + 
  scale_y_continuous("Cummulative number of behavior", breaks = seq(0, 60, by=2)) + 
  scale_x_discrete("Observation date")

# Leipzig B
Leipzig <- Data[,c("ObservationID","Group","Behavior")] 
LeipzigB <- Leipzig[(Leipzig$Group=="B"),] 
LeipzigB$Date <- str_sub(LeipzigB$ObservationID, start = 1, end = 4)
BbehaviorsID <- aggregate(data=LeipzigB, Behavior ~ Date, function(x) length(unique(x))) ; colnames(BbehaviorsID) <- c("Date","N")
LeipzigB <- merge(LeipzigB, BbehaviorsID, by=c("Date"))
BRepertoireID <- distinct(LeipzigB,Behavior,Date)

LeipzigB2 <- as.data.table(unique(BRepertoireID))
setkey(LeipzigB2, "Date")
LeipzigB2[, Behavior := as.numeric(factor(Behavior, levels = unique(Behavior)))]
setkey(LeipzigB2, "Date", "Behavior")
LeipzigB2.out <- LeipzigB2[J(unique(Date)), mult="last"]
LeipzigB2.out[, Behavior := cummax(Behavior)]

ggplot(LeipzigB2.out, aes(Date, Behavior)) + 
  geom_line(aes(group = 1)) +
  geom_point(size=3) + 
  theme_angele_ss + 
  scale_y_continuous("Cummulative number of behavior", breaks = seq(0, 60, by=2)) + 
  scale_x_discrete("Observation date")

# C2
C2 <- MediaFile[(MediaFile$Group=="C2"),] 
C2$`Date` <- dirname(C2$`Media file`)
write_xlsx(C2,"C:/Users/angy9/OneDrive - UT Cloud/Documents/STATS/Article 2/C2.xlsx")
#Delete "/Goup 2" at end of path + import C2
C2$`Date` <- str_sub(C2$`Date`, -4, -1)
C2RepertoireID <- distinct(C2,Behavior,Date)
C2_2 <- as.data.table(unique(C2RepertoireID))
setkey(C2_2, "Date")
C2_2[, Behavior := as.numeric(factor(Behavior, levels = unique(Behavior)))]
setkey(C2_2, "Date", "Behavior")
C2_2.out <- C2_2[J(unique(Date)), mult="last"]
C2_2.out[, Behavior := cummax(Behavior)]

ggplot(C2_2.out, aes(Date, Behavior)) + 
  geom_line(aes(group = 1)) +
  geom_point(size=3) + 
  theme_angele_ss + 
  scale_y_continuous("Cummulative number of behavior", breaks = seq(0, 60, by=2)) + 
  scale_x_discrete("Observation date")

# C4
C4 <- MediaFile[(MediaFile$Group=="C4"),]
C4$`Date` <- dirname(C4$`Media file`)
write_xlsx(C4,"C:/Users/angy9/OneDrive - UT Cloud/Documents/STATS/Article 2/C4.xlsx")
# Delete "/Goup 4" at end of path + import C4
C4$`Date` <- str_sub(C4$`Date`, -4, -1)
C4RepertoireID <- distinct(C4,Behavior,Date)
C4_2 <- as.data.table(unique(C4RepertoireID))
setkey(C4_2, "Date")
C4_2[, Behavior := as.numeric(factor(Behavior, levels = unique(Behavior)))]
setkey(C4_2, "Date", "Behavior")
C4_2.out <- C4_2[J(unique(Date)), mult="last"]
C4_2.out[, Behavior := cummax(Behavior)]

ggplot(C4_2.out, aes(Date, Behavior)) + 
  geom_line(aes(group = 1)) +
  geom_point(size=3) + 
  theme_angele_ss + 
  scale_y_continuous("Cummulative number of behavior", breaks = seq(0, 60, by=2)) + 
  scale_x_discrete("Observation date")

# Leintalzoo
Leintalzoo <- MediaFile[(MediaFile$Group=="L"),]
Leintalzoo$`Date` <- dirname(Leintalzoo$`Media file`)
Leintalzoo$`Date` <- str_sub(Leintalzoo$`Date`, -4, -1)
LeintalzooRepertoireID <- distinct(Leintalzoo,Behavior,Date)
Leintalzoo2 <- as.data.table(unique(LeintalzooRepertoireID))
setkey(Leintalzoo2, "Date")
Leintalzoo2[, Behavior := as.numeric(factor(Behavior, levels = unique(Behavior)))]
setkey(Leintalzoo2, "Date", "Behavior")
Leintalzoo2.out <- Leintalzoo2[J(unique(Date)), mult="last"]
Leintalzoo2.out[, Behavior := cummax(Behavior)]

ggplot(Leintalzoo2.out, aes(Date, Behavior)) + 
  geom_line(aes(group = 1)) +
  geom_point(size=3) + 
  theme_angele_ss + 
  scale_y_continuous("Cummulative number of behavior", breaks = seq(0, 60, by=2)) + 
  scale_x_discrete("Observation date")

# Kanyawara
Kanyawara <- MediaFile[(MediaFile$Group=="K"),]
Kanyawara$`Date` <- str_sub(Kanyawara$`Media file`, 31, -1)
Kanyawara$`Year` <- str_sub(Kanyawara$`Date`, 1, 4)
Kanyawara <- Kanyawara %>%
  mutate(Date = if_else(Year == 2013,
                        str_sub(Date, 6, 15),
                        str_sub(Date, 9, 18)))
Kanyawara <- Kanyawara %>%
  mutate(Date = if_else(
    Year == 2013,
    paste0(str_sub(Date, 1, 2), str_sub(Date, 4, 5), str_sub(Date, 7, nchar(Date))),
    paste0(str_sub(Date, 1, 4), str_sub(Date, 6, 7), str_sub(Date, 9, nchar(Date)))
  ))
Kanyawara <- Kanyawara %>%
  mutate(Date = if_else(
    Year == 2013,
    paste0(
      str_sub(Date, 5, 8),  # YYYY
      str_sub(Date, 3, 4),  # MM
      str_sub(Date, 1, 2)   # DD
    ),
    Date  # Leave other years unchanged
  ))
KanyawaraRepertoireID <- distinct(Kanyawara,Behavior,Date)
Kanyawara2 <- as.data.table(unique(KanyawaraRepertoireID))
setkey(Kanyawara2, "Date")
Kanyawara2[, Behavior := as.numeric(factor(Behavior, levels = unique(Behavior)))]
setkey(Kanyawara2, "Date", "Behavior")
Kanyawara2.out <- Kanyawara2[J(unique(Date)), mult="last"]
Kanyawara2.out[, Behavior := cummax(Behavior)]

ggplot(Kanyawara2.out, aes(Date, Behavior)) + 
  geom_line(aes(group = 1)) +
  geom_point(size=3) + 
  theme_angele_ss + 
  scale_y_continuous("Cummulative number of behavior", breaks = seq(0, 60, by=2)) + 
  scale_x_discrete("Observation date")

# Put them all on the same graph
LeipzigA2.out$Group <- "A" ; LeipzigA2.out<-LeipzigA2.out[,c("Behavior","Date","Group")] ; LeipzigA2.out$Date <- c(1:39)
LeipzigB2.out$Group <- "B" ; LeipzigB2.out<-LeipzigB2.out[,c("Behavior","Date","Group")] ; LeipzigB2.out$Date <- c(1:31)
C2_2.out$Group <- "C2" ; C2_2.out$Date <- c(1:49)
C4_2.out$Group <- "C4" ; C4_2.out$Date <- c(1:24)
Leintalzoo2.out$Group <- "L" ; Leintalzoo2.out$Date <- c(1:43)
Kanyawara2.out$Group <- "K" ; Kanyawara2.out$Date <- c(1:150)
CumulativeRepertoire <- rbind(LeipzigA2.out, LeipzigB2.out, C2_2.out, C4_2.out, Leintalzoo2.out, Kanyawara2.out)
CumulativeRepertoire$Group <- factor(CumulativeRepertoire$Group,levels = c("A", "B", "L", "C2", "C4", "K"))
ggplot(CumulativeRepertoire, aes(x=Date, y=Behavior, group=Group, color=Group)) + 
  geom_line(linewidth=1.5) +
  geom_point(aes(color=Group),size=2) + 
  theme_angele_ss + 
  scale_color_manual(values = c("#007ABB","#00AFBB", "#259C39", "#E7B800","#E79A00", "#ba0f09")) + 
  scale_y_continuous("Cummulative number of observed signals", breaks = seq(0, 80, by=10)) + 
  scale_x_continuous("Observation time (days)", breaks = seq(0, 150, by=10))

_____________________________________________________________________________________________________________________
2. Repertoire size 

# Individual repertoire size 
IndRepSize <- aggregate(data=Data, Behavior ~ Subject, function(x) length(unique(x))) ; colnames(IndRepSize) <- c("Subject","RepertoireSize") 

GroupRepSize <- aggregate(data=Data, Behavior ~ Group, function(x) length(unique(x))) ; colnames(GroupRepSize) <- c("Group","RepertoireSize") 


##Graph
IDinfo <- Data[,c("Subject","Group","Setting","Sex","Age","Rank", "GroupSize","NbAdultMales","NbJuveniles","NbInteractionPartners","TotComActs","ASOSize")] # extraire les colonnes du tableau d'origine
IDinfo <- IDinfo %>%
  group_by(across(-Age)) %>%
  slice_max(order_by = Age, n = 1, with_ties = FALSE) %>%
  ungroup()
write_xlsx(IDinfo,"C:/Users/angy9/OneDrive - UT Cloud/Documents/STATS/Article 2/IDinfo.xlsx")

RepComp <- merge(IndRepSize, IDinfo, by.x = "Subject", by.y = "Subject", all.x = TRUE, all.y = TRUE)
RepComp$Group <- factor(RepComp$Group,levels = c("A", "B", "L", "C2", "C4", "K"))

library(tidytext)
ggplot(RepComp, aes(x = reorder(Subject, RepertoireSize), y = RepertoireSize)) + 
  geom_col(linewidth = 1.5, aes(fill = Group, color = Sex)) + 
  theme_angele_ss + 
  scale_x_discrete("Subject") +
  scale_y_continuous("Repertoire size", breaks = seq(0, 30, by = 5)) + 
  scale_fill_manual(values = c("#007ABB", "#00AFBB","#259C39", "#E7B800","#E79A00", "#ba0f09")) + 
  scale_color_manual(values = c("white", "black")) +
  theme(legend.text = element_text(vjust = 2.75), axis.text.x = element_blank())  # OR axis.text.x = element_blank()

RepComp <- RepComp |> arrange(Group, RepertoireSize) |> mutate(Subject_ord = factor(Subject, levels = Subject))
ggplot(RepComp, aes(x = Subject_ord, y = RepertoireSize)) +
  geom_col(aes(fill = Group, color = Sex), linewidth = 1.5) +
  theme_angele_ss +
  scale_x_discrete("Subject") +
  scale_y_continuous("Repertoire size", breaks = seq(0, 30, by = 5)) +
  scale_fill_manual(values = c("#007ABB", "#00AFBB", "#259C39", "#E7B800", "#E79A00", "#ba0f09")) +
  scale_color_manual(values = c("white", "black")) +
  theme(legend.text = element_text(vjust = 2.75), axis.text.x = element_blank())


##MODEL 

RepComp$Age <- as.numeric(RepComp$Age)

vif(lm(log(RepertoireSize) ~ Age + Sex + Setting + NbInteractionPartners + NbAdultMales + log(TotComActs) + ASOSize, data = RepComp)) ### collinearity check ok: max vif = 1.7 < 5 
#                          GVIF Df GVIF^(1/(2*Df))
#Age                   1.142824  1        1.069029
#Sex                   1.090579  1        1.044308
#Setting               3.103037  2        1.327232
#NbInteractionPartners 2.613195  1        1.616538
#NbAdultMales          3.105752  1        1.762315
#log(TotComActs)       1.901486  1        1.378944
#ASOSize               1.945425  1        1.394785

plot.mod.RepertoireSize.m_null <- lmer(formula = log(RepertoireSize) ~ Age + Sex + log(TotComActs) + (1|Group), data = RepComp)

plot.mod.RepertoireSize.m <- lmer(formula = log(RepertoireSize) ~ Age + Sex + Setting + NbInteractionPartners + NbAdultMales + log(TotComActs) + ASOSize + (1|Group), data = RepComp)

as.data.frame(anova(plot.mod.RepertoireSize.m_null, plot.mod.RepertoireSize.m, test="LTR"))
#                               npar        AIC      BIC    logLik  -2*log(L)    Chisq Df   Pr(>Chisq)
#plot.mod.RepertoireSize.m_null    6   3.906997 19.10259  4.046501  -8.093003       NA NA           NA
#plot.mod.RepertoireSize.m        11 -11.838316 16.02028 16.919158 -33.838316 25.74531  5 9.997853e-05 --> p-value<0.05 so full model is better than null model

summary(plot.mod.RepertoireSize.m)
# Random effects:
#Groups   Name        Variance Std.Dev.
#Group    (Intercept) 0.03643  0.1909  
#Residual             0.04078  0.2019  
#Number of obs: 93, groups:  Group, 6

# Fixed effects:
#                       Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)            0.572029   0.280784 13.643389   2.037 0.061508 .  
#Age                   -0.002321   0.002397 81.964537  -0.968 0.335874    
#SexM                  -0.020670   0.051247 81.935033  -0.403 0.687740    
#SettingSanctuary      -0.086518   0.186989  1.679735  -0.463 0.696492    
#SettingWild           -0.277370   0.302189  1.515111  -0.918 0.480931    
#NbInteractionPartners  0.002768   0.005704 72.093035   0.485 0.628902    
#NbAdultMales           0.001619   0.021877  1.931590   0.074 0.947956    
#log(TotComActs)        0.372904   0.051917 83.481701   7.183 2.62e-10 ***
#ASOSize                0.046560   0.012483 83.782632   3.730 0.000347 ***

drop1(plot.mod.RepertoireSize.m,  test ="Chisq")
#                       Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Age                   0.03821 0.03821     1 81.965  0.9371 0.3358737    
#Sex                   0.00663 0.00663     1 81.935  0.1627 0.6877400    
#Setting               0.03579 0.01789     2  2.000  0.4388 0.6950142    
#NbInteractionPartners 0.00961 0.00961     1 72.093  0.2356 0.6289020    
#NbAdultMales          0.00022 0.00022     1  1.932  0.0055 0.9479561    
#log(TotComActs)       2.10386 2.10386     1 83.482 51.5917 2.618e-10 ***
#ASOSize               0.56733 0.56733     1 83.783 13.9123 0.0003475 ***
  
rep_plot.mod.RepertoireSize.m <- rpt(log(RepertoireSize) ~ Age + Sex + Setting + NbInteractionPartners + NbAdultMales + log(TotComActs) + ASOSize + (1|Group), 
              grname = "Group", data = RepComp, datatype = "Gaussian", nboot = 1000, npermut = 0, adjusted = FALSE)
print(rep_plot.mod.RepertoireSize.m)
#Repeatability estimation using the lmm method 

#Repeatability for Group
#R  = 0.194
#SE = 0.149
#CI = [0, 0.507]
#P  = 1 [LRT]
#NA [Permutation]


#PLOT
library(ggeffects)

predRepComp <- ggpredict(plot.mod.RepertoireSize.m, terms = "ASOSize")

RepPlot <- ggplot(predRepComp, aes(x = x, y = predicted)) +
  geom_point(data = RepComp, aes(x = ASOSize, y = RepertoireSize, size = TotComActs, col = factor(Group)), shape = 1, show.legend = FALSE) +
  geom_line(color = "black",  size=1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "#69b3a2") +
  theme_angele_ss +
  ggtitle("(a)") +
  scale_x_continuous("Number of social goals", breaks = seq(min(RepComp$ASOSize), max(RepComp$ASOSize), by = 1)) +
  scale_y_continuous("Repertoire size", breaks = seq(0, 30, by = 5)) +
  scale_color_manual(name="Group: ", values = c("#007ABB", "#00AFBB", "#259C39", "#E7B800", "#E79A00", "#ba0f09")) +
  theme(plot.margin = margin(t = 4, b = 6, r = 10, l = 4)) 




