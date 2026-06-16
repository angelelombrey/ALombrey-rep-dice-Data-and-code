# PACKAGES
library(lme4);library(tidyverse); library(dplyr);library(writexl);library(ggpubr);library(lmerTest)
library(plyr);library(tidyr);library(car);library(RColorBrewer);library(ggplot2)

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
                         strip.text = element_text(size = 15))

_______________________________________________________________________________________________________________________________________________________________________________________________
# PREPARE DATA

Data$SexRec <- str_sub(Data$SexCombinaison, start = 2)
Data <- Data %>% unite("SexSetting", c(Sex,Setting), sep = "/", remove = FALSE)

# Take out the Keepers-directed signals, None and NA in Recipient
Data <- Data[Data$Recipient != "Keepers", ]; Data <- Data[Data$Recipient != "NA", ]; Data <- Data[Data$Recipient != "None", ]
Data <- subset(Data, !is.na(Recipient))
# Take out lines "no focal subject"
Data <- Data[Data$Subject != "No focal subject", ] 

Data <- Data %>% mutate(across(c(Behavior2, Behavior3, Behavior4, Behavior5), str_trim))
Data$Behavior2 <- str_to_sentence(Data$Behavior2); Data$Behavior3 <- str_to_sentence(Data$Behavior3); Data$Behavior4 <- str_to_sentence(Data$Behavior4); Data$Behavior5 <- str_to_sentence(Data$Behavior5)

# Take out "Special" in Behaviors
Data <- Data[Data$Behavior != "Special", ] ; Data <- Data[is.na(Data$Behavior2) | Data$Behavior2 != "Special", ]
Data$Behavior[Data$Behavior == "Pout"] <- "Pant hoot face" ; Data$Behavior2[Data$Behavior2 == "Pout"] <- "Pant hoot face" ; Data$Behavior3[Data$Behavior3 == "Pout"] <- "Pant hoot face" ; Data$Behavior4[Data$Behavior4 == "Pout"] <- "Pant hoot face" ; Data$Behavior5[Data$Behavior5 == "Pout"] <- "Pant hoot face"
Data$Behavior[Data$Behavior == "Present grooming"] <- "Present" ; Data$Behavior2[Data$Behavior2 == "Present grooming"] <- "Present" ; Data$Behavior3[Data$Behavior3 == "Present grooming"] <- "Present" ; Data$Behavior4[Data$Behavior4 == "Present grooming"] <- "Present" ; Data$Behavior5[Data$Behavior5 == "Present grooming"] <- "Present"
Data$Behavior[Data$Behavior == "Present sexual"] <- "Present" ; Data$Behavior2[Data$Behavior2 == "Present sexual"] <- "Present" ; Data$Behavior3[Data$Behavior3 == "Present sexual"] <- "Present" ; Data$Behavior4[Data$Behavior4 == "Present sexual"] <- "Present" ; Data$Behavior5[Data$Behavior5 == "Present sexual"] <- "Present" 
Data$Behavior[Data$Behavior == "Present climb on me"] <- "Present" ; Data$Behavior2[Data$Behavior2 == "Present climb on me"] <- "Present" ; Data$Behavior3[Data$Behavior3 == "Present climb on me"] <- "Present" ; Data$Behavior4[Data$Behavior4 == "Present climb on me"] <- "Present" ; Data$Behavior5[Data$Behavior5 == "Present climb on me"] <- "Present"

library(dplyr)
library(purrr)
library(stringr)

# Multicomponent Combinations
Data <- Data %>%
  mutate(
    SignalCombination = pmap_chr(
      list(Behavior, Behavior2, Behavior3, Behavior4, Behavior5),
      function(b1, b2, b3, b4, b5) {
        vals <- c(b1, b2, b3, b4, b5) %>%
          discard(~ is.na(.) || . == "")
        
        # if only Behavior is present → NA
        if (length(vals) <= 1) {
          return(NA_character_)
        }
        
        sort(vals) %>% str_c(collapse = " + ")
      }
    )
  )

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

# Take out duplicated lines
Datareduit <- Data[Data$Duplicated != "T", ]  
datMultiRep <- Datareduit

# Sampling effort
SamplingEffort <- data.frame(datMultiRep %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(SamplingEffort) <- c("Subject","TotComActs") 
datMultiRep <- merge(SamplingEffort, datMultiRep, by=c("Subject"))

datMultiRep <- datMultiRep[!is.na(datMultiRep$Group), ]

datMultiRep <- datMultiRep[datMultiRep$Multimodality_YN != "0", ]; datMultiRep <- datMultiRep[datMultiRep$Multimodality_YN != "NA", ]; datMultiRep <- subset(datMultiRep, !is.na(Multimodality_YN))
datMultiRep <- datMultiRep[!is.na(datMultiRep$SignalCombination), ]

# Individual ASO size 
datMultiRep <- datMultiRep[datMultiRep$ASO != "NA", ]; datMultiRep <- datMultiRep[!is.na(datMultiRep$ASO), ]
IndASOSize <- aggregate(data=datMultiRep, ASO ~ Subject, function(x) length(unique(x))) ; colnames(IndASOSize) <- c("Subject","ASOSize") 
datMultiRep <- merge(IndASOSize, datMultiRep, by=c("Subject"))

# Omit all levels of Subject that contributed fewer than 5 cases --- 118 --> 81 individuals (54F, 27M, 30C, 34S, 17W)
nb.signals.per.ind <- data.frame(datMultiRep %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(nb.signals.per.ind) <- c("Subject","TotComActs") 
nb.signals.per.ind$Subject <- as.character(nb.signals.per.ind$Subject)
nb.signals.per.ind <- subset(nb.signals.per.ind, nb.signals.per.ind$TotComActs>4) # Keeps only the values in "column" that appear more than 4 times (N >= 4)
datMultiRep <- subset(datMultiRep, Subject %in% nb.signals.per.ind$Subject)


# Multi Sampling effort
MultiSamplingEffort <- data.frame(datMultiRep %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(MultiSamplingEffort) <- c("Subject","TotMultiSignals") 
datMultiRep <- merge(MultiSamplingEffort, datMultiRep, by=c("Subject"))


length(unique(datMultiRep[["SignalCombination"]])) #406

________________________________________________________________________________________________________________________________________________________________________________________________
3. Multimodal Repertoire size 

# Individual multimodal repertoire size 
IndMultiRepSize <- aggregate(data=datMultiRep, SignalCombination ~ Subject, function(x) length(unique(x))) ; colnames(IndMultiRepSize) <- c("Subject","MultiRepSize") 

##Graph
IDinfo <- datMultiRep[,c("Subject","Group","Setting","Sex","Age","Rank", "GroupSize","NbAdultMales","NbJuveniles","NbInteractionPartners","TotComActs","TotMultiSignals","ASOSize")] 
IDinfo <- IDinfo %>%
  group_by(across(-Age)) %>%
  slice_max(order_by = Age, n = 1, with_ties = FALSE) %>%
  ungroup()
MultiRepComp <- merge(IndMultiRepSize, IDinfo, by.x = "Subject", by.y = "Subject", all.x = TRUE, all.y = TRUE)
MultiRepComp$Group <- factor(MultiRepComp$Group,levels = c("A", "B", "L", "C2", "C4", "K"))

ggplot(MultiRepComp, aes(x = reorder(Subject, MultiRepSize), y = MultiRepSize)) + 
  geom_col(linewidth = 1.5, aes(fill = Group, color = Sex)) + 
  theme_angele_ss + 
  scale_x_discrete("Subject") +
  scale_y_continuous("Multimodal Repertoire size", breaks = seq(0, 40, by=5), limits = c(0, 40)) + 
  scale_fill_manual(values = c("#007ABB", "#00AFBB","#259C39", "#E7B800","#E79A00", "#ba0f09")) + 
  scale_color_manual(values = c("white", "black")) +
  theme(legend.text = element_text(vjust = 2), axis.text.x = element_blank())  # OR axis.text.x = element_blank()


##MODEL

MultiRepComp$Age <- as.numeric(MultiRepComp$Age)

vif(lm(log(MultiRepSize) ~ Age + Sex + Setting + NbInteractionPartners + NbAdultMales + log(TotComActs) + ASOSize, data = MultiRepComp)) ### collinearity check ok: max vif = 1.7 < 5 
#                          GVIF Df GVIF^(1/(2*Df))
#Age                   1.034246  1        1.016979
#Sex                   1.206030  1        1.098194
#Setting               2.919946  2        1.307205
#NbInteractionPartners 2.484484  1        1.576225
#NbAdultMales          3.083759  1        1.756063
#log(TotComActs)       1.892492  1        1.375679
#ASOSize               1.314844  1        1.146667

plot.mod.MultiRepSize.m_null <- lmer(formula = log(MultiRepSize) ~ Age + Sex + log(TotComActs) + (1|Group), data = MultiRepComp)

plot.mod.MultiRepSize.m <- lmer(formula = log(MultiRepSize) ~ Age + Sex + Setting + NbInteractionPartners + NbAdultMales + log(TotComActs) + ASOSize + (1|Group), data = MultiRepComp)

as.data.frame(anova(plot.mod.MultiRepSize.m_null, plot.mod.MultiRepSize.m, test="LTR"))
#                             npar      AIC       BIC    logLik  -2*log(L)    Chisq Df   Pr(>Chisq)
#plot.mod.MultiRepSize.m_null    6 94.31354 108.68023 -41.15677  82.31354        NA NA           NA
#plot.mod.MultiRepSize.m        11 57.16154  83.50048 -17.58077  35.16154    47.152  5 5.290178e-09 --> p-value<0.05 so full model is better than null model

summary(plot.mod.MultiRepSize.m)
# Random effects:
#Groups   Name        Variance Std.Dev.
#Group    (Intercept) 0.04997  0.2235  
#Residual             0.09511  0.3084  
#Number of obs: 81, groups:  Group, 6

# Fixed effects:
#                       Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)           -0.057086   0.401279 21.033831  -0.142    0.888    
#Age                    0.002507   0.004008 71.506103   0.626    0.534    
#SexM                  -0.147581   0.083105 71.996936  -1.776    0.080 .  
#SettingSanctuary      -0.017682   0.231957  1.565683  -0.076    0.948    
#SettingWild           -0.266541   0.364292  1.300394  -0.732    0.573    
#NbInteractionPartners -0.001298   0.008263 63.090644  -0.157    0.876    
#NbAdultMales          -0.009819   0.027931  1.984964  -0.352    0.759    
#log(TotComActs)        0.421855   0.082377 71.167418   5.121 2.49e-06 ***
#ASOSize                0.135560   0.019960 70.144306   6.791 2.97e-09 ***

drop1(plot.mod.MultiRepSize.m,  test ="Chisq")
#                      Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Age                   0.0372  0.0372     1 71.506  0.3913   0.53361    
#Sex                   0.3000  0.3000     1 71.997  3.1536   0.07999 .  
#Setting               0.0528  0.0264     2  2.000  0.2776   0.78273    
#NbInteractionPartners 0.0023  0.0023     1 63.091  0.0247   0.87565    
#NbAdultMales          0.0118  0.0118     1  1.985  0.1236   0.75897    
#log(TotComActs)       2.4944  2.4944     1 71.167 26.2247 2.491e-06 ***
#ASOSize               4.3871  4.3871     1 70.144 46.1244 2.966e-09 ***
  
rep_plot.mod.MultiRepSize.m <- rpt(log(MultiRepSize) ~ Age + Sex + Setting + NbInteractionPartners + NbAdultMales + log(TotComActs) + ASOSize + (1|Group), 
                                     grname = "Group", data = MultiRepComp, datatype = "Gaussian", nboot = 1000, npermut = 0, adjusted = FALSE)
print(rep_plot.mod.MultiRepSize.m)
#Repeatability estimation using the lmm method 

#Repeatability for Group
#R  = 0.128
#SE = 0.12
#CI = [0, 0.401]
#P  = 0.5 [LRT]
# NA [Permutation]


# PLOT
predMultiRepComp <- ggpredict(plot.mod.MultiRepSize.m, terms = "ASOSize")

MultiRepPlot <- ggplot(predMultiRepComp, aes(x = x, y = predicted)) +
  geom_point(data = MultiRepComp, aes(x = ASOSize, y = MultiRepSize, size = TotComActs, col = factor(Group)), shape = 1) +
  geom_line(color = "black",  size=1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "#69b3a2") +
  theme_angele_ss +
  ggtitle("(b)") +
  scale_x_continuous("Number of social goals", breaks = seq(min(MultiRepComp$ASOSize), max(RepComp$ASOSize), by = 1)) +
  scale_y_continuous("Multicomponent repertoire size", breaks = seq(0, 40, by = 5)) +
  scale_color_manual(name="Group: ", values = c("#007ABB", "#00AFBB", "#259C39", "#E7B800", "#E79A00", "#ba0f09")) +
  theme(legend.position = "right", legend.text = element_text(vjust = 4), plot.margin = margin(t = 4, b = 6, l = 10)) +
  guides(size = "none", color = guide_legend(override.aes = list(size = 4, alpha = 1)))

ggarrange(RepPlot, MultiRepPlot,widths = c(1, 1.2))




