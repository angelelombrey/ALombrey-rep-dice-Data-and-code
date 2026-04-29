5. SETTING MULTICOMPONENT REPERTOIRE SIMILARITY 

_______________________________________________________________________________________________________________________________________________________________________________________________
# PACKAGES
library(lme4);library(tidyverse); library(dplyr);library(writexl);library(ggpubr)
library(plyr);library(tidyr);library(car);library(RColorBrewer);library(ggplot2)
library(RColorBrewer)

rm(list=ls())

dodge.posn <- position_dodge(.9)
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
                         title = element_text(size = 25, family="sans"),
                         strip.text = element_text(size = 15))
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


# Take out duplicated lines
Datareduit <- Data[Data$Duplicated != "T", ]  
datMultiRep <- Datareduit

# Sampling effort
SamplingEffort <- data.frame(datMultiRep %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(SamplingEffort) <- c("Subject","TotComActs") 
datMultiRep <- merge(SamplingEffort, datMultiRep, by=c("Subject"))

datMultiRep <- datMultiRep[!is.na(datMultiRep$Group), ]

datMultiRep <- datMultiRep[datMultiRep$Multimodality_YN != "0", ]; datMultiRep <- datMultiRep[datMultiRep$Multimodality_YN != "NA", ]; datMultiRep <- subset(datMultiRep, !is.na(Multimodality_YN))
datMultiRep <- datMultiRep[!is.na(datMultiRep$SignalCombination), ]

# Omit all levels of Subject that contributed fewer than 5 cases --- 118 --> 19 individuals (16F, 3M, 0W)
nb.signals.per.ind <- data.frame(datMultiRep %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(nb.signals.per.ind) <- c("Subject","TotComActs") 
nb.signals.per.ind$Subject <- as.character(nb.signals.per.ind$Subject)
nb.signals.per.ind <- subset(nb.signals.per.ind, nb.signals.per.ind$TotComActs>4) # Keeps only the values in "column" that appear more than 4 times (N > 4)
datMultiRep <- subset(datMultiRep, Subject %in% nb.signals.per.ind$Subject)


# Multi Sampling effort
MultiSamplingEffort <- data.frame(datMultiRep %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(MultiSamplingEffort) <- c("Subject","TotMultiSignals") 
datMultiRep <- merge(MultiSamplingEffort, datMultiRep, by=c("Subject"))

IndMultiRepSize <- aggregate(data=datMultiRep, SignalCombination ~ Subject, function(x) length(unique(x))) ; colnames(IndMultiRepSize) <- c("Subject","MultiRepSize")
dat.ind.bhv.Multi2 <- data.frame(datMultiRep %>% group_by(Subject,SignalCombination) %>% dplyr::summarize(count=n())) ; colnames(dat.ind.bhv.Multi2) <- c("Subject","SignalCombination","NsignalsPbehavior") 

________________________________________________________________________________________________________________________________________________________________________________________________


# Individual repertoire similarity (Captive-Captive/Sanctuary-Sanctuary, Captive-Sanctuary): 

# Step 1 : calculate DiceMulti coefficient based on formula: dc = (2 x number of behaviours two inds have in common)/(R of ind1 + R ind2) ---
## create data frame in which every combination of individual and behaviour is counted
dat.ind.bhv.Multi <- data.frame(table(datMultiRep$Subject,datMultiRep$SignalCombination)) ; colnames(dat.ind.bhv.Multi) <- c("Subject","SignalCombination","Nsignals")
## create data frame with Subject as rows and columns
dat.ind.ind.Multi <- data.frame(matrix(rep(NA, nrow(IndMultiRepSize) * nrow(IndMultiRepSize)), ncol=nrow(IndMultiRepSize)))
colnames(dat.ind.ind.Multi) <- IndMultiRepSize$Subject
rownames(dat.ind.ind.Multi) <- IndMultiRepSize$Subject
## loop over all pairs of individuals
for(i in 1:(nrow(IndMultiRepSize))) {
  for(j in 1:(nrow(IndMultiRepSize))) {
    Ind1 <- subset(dat.ind.bhv.Multi2, Subject==IndMultiRepSize$Subject[i])
    Ind2 <- subset(dat.ind.bhv.Multi2, Subject==IndMultiRepSize$Subject[j])
    ovrlp <- length(which(Ind1$SignalCombination %in% Ind2$SignalCombination)) # number of behaviors 2 individuals have in common
    RInd1 <- subset(IndMultiRepSize, Subject==IndMultiRepSize$Subject[i])$MultiRepSize	# which is the same as nrow(Ind1)
    RInd2 <- subset(IndMultiRepSize, Subject==IndMultiRepSize$Subject[j])$MultiRepSize	# which is the same as nrow(Ind2)
    dat.ind.ind.Multi[which(rownames(dat.ind.ind.Multi)==IndMultiRepSize$Subject[i]),which(colnames(dat.ind.ind.Multi)==IndMultiRepSize$Subject[j])] <- (2 * ovrlp) / (RInd1 + RInd2)
  }
}
## set diagonal to NA
DiceMulti <- as.matrix(dat.ind.ind.Multi)
diag(DiceMulti) <- NA
#DiceMulti <- as.data.frame(DiceMulti)

# step 2: define within- and between-setting dyads 
## create data frame with Subject as rows and columns
SETdat.wth.btw.Multi <- data.frame(matrix(rep(NA, nrow(IndMultiRepSize) * nrow(IndMultiRepSize)), ncol=nrow(IndMultiRepSize)))
colnames(SETdat.wth.btw.Multi) <- IndMultiRepSize$Subject
rownames(SETdat.wth.btw.Multi) <- IndMultiRepSize$Subject
## loop over all pairs of individuals
for(i in 1:(nrow(IndMultiRepSize))) {
  for(j in 1:(nrow(IndMultiRepSize))) {
    Ind1 <- subset(datMultiRep, Subject==IndMultiRepSize$Subject[i])
    Ind2 <- subset(datMultiRep, Subject==IndMultiRepSize$Subject[j])
    if(length(table(Ind1$Setting))>1) { print("Ohoh, better check my data again!") }
    if(length(table(Ind2$Setting))>1) { print("Ohoh, better check my data again!") }
    if(Ind1$Setting[1]==Ind2$Setting[1]) { SETwth.btw.Multi <- "Within" } else { SETwth.btw.Multi <- "Between" }
    SETdat.wth.btw.Multi[which(rownames(SETdat.wth.btw.Multi)==IndMultiRepSize$Subject[i]),which(colnames(SETdat.wth.btw.Multi)==IndMultiRepSize$Subject[j])] <- SETwth.btw.Multi
  }
}
## set diagonal to NA
SETwth.btw.Multi <- as.matrix(SETdat.wth.btw.Multi)
diag(SETwth.btw.Multi) <- NA
#SETwth.btw.Multi <- as.data.frame(SETwth.btw.Multi)

# step 3: select within and between setting dyads from the DiceMulti data frame (N=19)
DiceMulti.WithinSET <- ifelse(SETwth.btw.Multi=="Within",DiceMulti,NA)
DiceMulti.BetweenSET <- ifelse(SETwth.btw.Multi=="Between",DiceMulti,NA)

(Emp.MDiceMultiWithinSET <- mean(DiceMulti.WithinSET[lower.tri(DiceMulti.WithinSET, diag = FALSE)],na.rm=TRUE)) # mean DiceMulti coefficient for within-dyCMultids (= mean repertoire overlap between individuals in the same setting (within, A-A and B-B))
# Mean DiceMulti within settings (settings C, S et W confondus) = 0.07943728
(Emp.MDiceMultiBetweenSET <- mean(DiceMulti.BetweenSET[lower.tri(DiceMulti.BetweenSET, diag = FALSE)],na.rm=TRUE)) # mean DiceMulti coefficient for between-dyCMultids (= mean repertoire overlap between individuals in different setting (between, A-B))
# Mean DiceMulti between (C-S) settings = 0.04682856


# Repertoire similarity between groups Within et Between

# step 4: Combine all info in 1 table (1 line = 1 dyad --> create data frame dd with columns Ind1, Ind2, Dice, Within/Between)
## Tranform Within matrix in table
yWithinSETMulti <- expand.grid(rownames(DiceMulti.WithinSET), colnames(DiceMulti.WithinSET)) 
labsWithin <- yWithinSETMulti[as.vector(upper.tri(DiceMulti.WithinSET, diag = FALSE)), ] 
yWithinSETMulti <- cbind(labsWithin, DiceMulti.WithinSET[upper.tri(DiceMulti.WithinSET,diag=FALSE)]) 
colnames(yWithinSETMulti) <- c("Ind1","Ind2","DiceMulti")
yWithinSETMulti <- yWithinSETMulti[!is.na(yWithinSETMulti$DiceMulti), ] 
yWithinSETMulti$WithinBetween <- rep("WithinS",nrow(yWithinSETMulti)) 
mean(yWithinSETMulti$Dice)

## Tranform Between matrix in table
yBetweenSETMulti <- expand.grid(rownames(DiceMulti.BetweenSET), colnames(DiceMulti.BetweenSET))
labsBetween <- yBetweenSETMulti[as.vector(upper.tri(DiceMulti.BetweenSET, diag = FALSE)), ]
yBetweenSETMulti <- cbind(labsBetween, DiceMulti.BetweenSET[upper.tri(DiceMulti.BetweenSET,diag=FALSE)])
colnames(yBetweenSETMulti) <- c("Ind1","Ind2","DiceMulti")
yBetweenSETMulti <- yBetweenSETMulti[!is.na(yBetweenSETMulti$DiceMulti), ]
yBetweenSETMulti$WithinBetween <- rep("BetweenS",nrow(yBetweenSETMulti))
mean(yBetweenSETMulti$Dice)

##Combine Within and Between tables
WBddMultiSET <- rbind(yWithinSETMulti,yBetweenSETMulti)

datMultiRep.sub <- datMultiRep[,c("Subject","Setting")] 
datMultiRep.sub <- datMultiRep.sub[!duplicated(datMultiRep.sub), ] 
WBddMultiSET <- merge(WBddMultiSET, datMultiRep.sub, by.x="Ind1", by.y="Subject", sort=FALSE, all.y=FALSE) 
WBddMultiSET <- merge(WBddMultiSET, datMultiRep.sub, by.x="Ind2", by.y="Subject", sort=FALSE) 
colnames(WBddMultiSET)[which(colnames(WBddMultiSET)=="Setting.x")] <- "Setting.ind1"
colnames(WBddMultiSET)[which(colnames(WBddMultiSET)=="Setting.y")] <- "Setting.ind2"
# Order columns
WBddMultiSET<-WBddMultiSET[,c("Ind1","Ind2","DiceMulti","WithinBetween","Setting.ind1","Setting.ind2")]

# step 5: matrix permutation test ---
nPerm <- 1000
dat.Perm.MultiSET <- c()
for(k in 1:nPerm) {
  ## permute Setting in datMultiRep
  ind.Setting.perm.Multi <- datMultiRep.sub
  ind.Setting.perm.Multi$Setting <- sample(ind.Setting.perm.Multi$Setting, nrow(ind.Setting.perm.Multi), replace=FALSE)
  datMultiRep.inc.perm.MultiSET <- datMultiRep[,-which(colnames(datMultiRep)=="Setting")] 
  datMultiRep.inc.perm.MultiSET <- merge(datMultiRep.inc.perm.MultiSET, ind.Setting.perm.Multi, by.x="Subject", by.y="Subject", sort=FALSE)
  ## repeat step 2 & 3
  ## create data frame with animal_id as rows and columns (step2)
  SETdat.wth.btw.Multi <- data.frame(matrix(rep(NA, nrow(IndMultiRepSize) * nrow(IndMultiRepSize)), ncol=nrow(IndMultiRepSize)))
  colnames(SETdat.wth.btw.Multi) <- IndMultiRepSize$Subject
  rownames(SETdat.wth.btw.Multi) <- IndMultiRepSize$Subject
  ## loop over all pairs of individuals (step2)
  for(i in 1:(nrow(IndMultiRepSize))) {
    for(j in 1:(nrow(IndMultiRepSize))) {
      Ind1 <- subset(datMultiRep.inc.perm.MultiSET, Subject==IndMultiRepSize$Subject[i])
      Ind2 <- subset(datMultiRep.inc.perm.MultiSET, Subject==IndMultiRepSize$Subject[j])
      if(length(table(Ind1$Setting))>1) { print("Ohoh, better check my data again!") }
      if(length(table(Ind2$Setting))>1) { print("Ohoh, better check my data again!") }
      if(Ind1$Setting[1]==Ind2$Setting[1]) { SETwth.btw.Multi <- "Within" } else { SETwth.btw.Multi <- "Between" }
      SETdat.wth.btw.Multi[which(rownames(SETdat.wth.btw.Multi)==IndMultiRepSize$Subject[i]),which(colnames(SETdat.wth.btw.Multi)==IndMultiRepSize$Subject[j])] <- SETwth.btw.Multi
    }
  }
  ## set diagonal to NA (step3)
  SETwth.btw.Multi <- as.matrix(SETdat.wth.btw.Multi)
  diag(SETwth.btw.Multi) <- NA
  DiceMulti.WithinSET <- ifelse(SETwth.btw.Multi=="Within",DiceMulti,NA)
  DiceMulti.BetweenSET <- ifelse(SETwth.btw.Multi=="Between",DiceMulti,NA)
  MDiceMultiWithin <- mean(DiceMulti.WithinSET[lower.tri(DiceMulti.WithinSET, diag = FALSE)],na.rm=TRUE)
  MDiceMultiBetween <- mean(DiceMulti.BetweenSET[lower.tri(DiceMulti.BetweenSET, diag = FALSE)],na.rm=TRUE)
  dat.Perm.MultiSET[k] <- MDiceMultiWithin - MDiceMultiBetween
  flush.console()
  if(k %% 10 == 0) { print(paste0("Finished ", k, " out of ", nPerm, " simulations")) } 
}
hist(dat.Perm.MultiSET)
abline(v=Emp.MDiceMultiWithinSET - Emp.MDiceMultiBetweenSET, col="red") 
pMRS <- (sum(abs(dat.Perm) >= abs(Emp.MDiceMultiWithin - Emp.MDiceMultiBetween)) + 1) / (nPerm + 1)
# P-value = 0.000999


# Significance thresholds are P>=0.975 and P<=0.025. 
# This is because we are looking at the deviation from 0, which can either be negative or positive.
# Because the distribution of differences is not necessarily symmetric around zero, we cannot calculate P-values using absolute values.


# Repertoire similarity within setting captive (C-C), within setting sanctuary (S-S) & within setting wild (W-W)

# step 2b: define within and between setting diads --- for C and S settings 
## create data frame with animal_id as rows and columns
SETdat.wth.btw.Multi2 <- data.frame(matrix(rep(NA, nrow(IndMultiRepSize) * nrow(IndMultiRepSize)), ncol=nrow(IndMultiRepSize)))
colnames(SETdat.wth.btw.Multi2) <- IndMultiRepSize$Subject
rownames(SETdat.wth.btw.Multi2) <- IndMultiRepSize$Subject
## loop over all pairs of individuals
for(i in 1:(nrow(IndMultiRepSize))) {
  for(j in 1:(nrow(IndMultiRepSize))) {
    SETwth.btw.Multi2 <- NA
    Ind1 <- subset(datMultiRep, Subject==IndMultiRepSize$Subject[i])
    Ind2 <- subset(datMultiRep, Subject==IndMultiRepSize$Subject[j])
    if(length(table(Ind1$Setting))>1) { print("Ohoh, better check my data again!") }
    if(length(table(Ind2$Setting))>1) { print("Ohoh, better check my data again!") }
    if(Ind1$Setting[1]=="Captive" & Ind2$Setting[1]=="Captive") { SETwth.btw.Multi2 <- "Captive" }
    if(Ind1$Setting[1]=="Sanctuary" & Ind2$Setting[1]=="Sanctuary") { SETwth.btw.Multi2 <- "Sanctuary" }
    if(Ind1$Setting[1]=="Wild" & Ind2$Setting[1]=="Wild") { SETwth.btw.Multi2 <- "Wild" }
    SETdat.wth.btw.Multi2[which(rownames(SETdat.wth.btw.Multi2)==IndMultiRepSize$Subject[i]),which(colnames(SETdat.wth.btw.Multi2)==IndMultiRepSize$Subject[j])] <- SETwth.btw.Multi2
  }
}
## set diagonal to NA
SETwth.btw.Multi2 <- as.matrix(SETdat.wth.btw.Multi2)
diag(SETwth.btw.Multi2) <- NA
#SETwth.btw.Multi2 <- as.data.frame(SETwth.btw.Multi2)

# step 3b: select within and between group diads from the DiceMulti data frame ---
DiceMulti.C <- ifelse(SETwth.btw.Multi2=="Captive",DiceMulti,NA)
DiceMulti.S <- ifelse(SETwth.btw.Multi2=="Sanctuary",DiceMulti,NA)
DiceMulti.W <- ifelse(SETwth.btw.Multi2=="Wild",DiceMulti,NA)

(Emp.MDiceMultiC <- mean(DiceMulti.C[lower.tri(DiceMulti.C, diag = FALSE)],na.rm=TRUE))
# Mean DiceMulti within Captive = 0.06299803 --> big individual variation (0 =  no overlap ; 1 = 100% overlap)
(Emp.MDiceMultiS <- mean(DiceMulti.S[lower.tri(DiceMulti.S, diag = FALSE)],na.rm=TRUE))
# Mean DiceMulti within Sanctuary = 0.08585487
(Emp.MDiceMultiW <- mean(DiceMulti.W[lower.tri(DiceMulti.W, diag = FALSE)],na.rm=TRUE))
# Mean DiceMulti within Wild = 0.1055461

## step 4b
## Tranform C matrix in table
yCMulti <- expand.grid(rownames(DiceMulti.C), colnames(DiceMulti.C)) 
labsC <- yCMulti[as.vector(upper.tri(DiceMulti.C, diag = FALSE)), ] 
yCMulti <- cbind(labsC, DiceMulti.C[upper.tri(DiceMulti.C,diag=FALSE)]) 
colnames(yCMulti) <- c("Ind1","Ind2","DiceMulti")
yCMulti <- yCMulti[!is.na(yCMulti$DiceMulti), ] 
yCMulti$SET <- rep("Captive",nrow(yCMulti))
mean(yCMulti$Dice)

## Tranform S matrix in table
ySMulti <- expand.grid(rownames(DiceMulti.S), colnames(DiceMulti.S))
labsS <- ySMulti[as.vector(upper.tri(DiceMulti.S, diag = FALSE)), ]
ySMulti <- cbind(labsS, DiceMulti.S[upper.tri(DiceMulti.S,diag=FALSE)])
colnames(ySMulti) <- c("Ind1","Ind2","DiceMulti")
ySMulti <- ySMulti[!is.na(ySMulti$DiceMulti), ]
ySMulti$SET <- rep("Sanctuary",nrow(ySMulti))
mean(ySMulti$Dice)

## Tranform W matrix in table
yWMulti <- expand.grid(rownames(DiceMulti.W), colnames(DiceMulti.W))
labsW <- yWMulti[as.vector(upper.tri(DiceMulti.W, diag = FALSE)), ]
yWMulti <- cbind(labsW, DiceMulti.W[upper.tri(DiceMulti.W,diag=FALSE)])
colnames(yWMulti) <- c("Ind1","Ind2","DiceMulti")
yWMulti <- yWMulti[!is.na(yWMulti$DiceMulti), ]
yWMulti$SET <- rep("Wild",nrow(yWMulti))
mean(yWMulti$Dice)

## Combine C, S and W tables
SETddMulti <- rbind(yCMulti,ySMulti,yWMulti)

datMultiRep.sub <- datMultiRep[,c("Subject","Setting")] 
datMultiRep.sub <- datMultiRep.sub[!duplicated(datMultiRep.sub), ] 
SETddMulti <- merge(SETddMulti, datMultiRep.sub, by.x="Ind1", by.y="Subject", sort=FALSE, all.y=FALSE) 
SETddMulti <- merge(SETddMulti, datMultiRep.sub, by.x="Ind2", by.y="Subject", sort=FALSE) 
colnames(SETddMulti)[which(colnames(SETddMulti)=="Setting.x")] <- "Setting.ind1"
colnames(SETddMulti)[which(colnames(SETddMulti)=="Setting.y")] <- "Setting.ind2"
# Order columns
SETddMulti<-SETddMulti[,c("Ind1","Ind2","DiceMulti","SET","Setting.ind1","Setting.ind2")]

mu.empMRS <- c(Emp.MDiceMultiC, Emp.MDiceMultiS, Emp.MDiceMultiW)
names(mu.empMRS) <- c("Captive","Sanctuary","Wild")
# Get the number of individuals per group
group_sizesSET <- table(datMultiRep.sub$Setting)[c("Captive","Sanctuary","Wild")]
wSET <- group_sizesSET * (group_sizesSET - 1) / 2
# Compute the weighted sums of squares
mean.w.empMRS <- weighted.mean(mu.empMRS, wSET=wSET, na.rm=TRUE)
SS.empMRS <- sum(wSET * (mu.empMRS - mean.w.empMRS)^2, na.rm=TRUE) 


# step 5b: matrix permutation test ---
nPerm <- 1000
dat.Perm.MultiSET2 <- c()
perm.pairwiseMRS <- matrix(NA_real_, nrow = nPerm, ncol = 3)
for(k in 1:nPerm) {
  ## permute group in datMultiRep
  ind.Setting.perm.Multi <- datMultiRep.sub
  ind.Setting.perm.Multi$Setting <- sample(ind.Setting.perm.Multi$Setting, nrow(ind.Setting.perm.Multi), replace=FALSE)
  datMultiRep.inc.perm.MultiSET <- datMultiRep[,-which(colnames(datMultiRep)=="Setting")] 
  datMultiRep.inc.perm.MultiSET <- merge(datMultiRep.inc.perm.MultiSET, ind.Setting.perm.Multi, by.x="Subject", by.y="Subject", sort=FALSE)
  ## repeat step 2 & 3
  ## create data frame with animal_id as rows and columns (step2)
  SETdat.wth.btw.Multi2 <- data.frame(matrix(rep(NA, nrow(IndMultiRepSize) * nrow(IndMultiRepSize)), ncol=nrow(IndMultiRepSize)))
  colnames(SETdat.wth.btw.Multi2) <- IndMultiRepSize$Subject
  rownames(SETdat.wth.btw.Multi2) <- IndMultiRepSize$Subject
  ## loop over all pairs of individuals (step2)
  for(i in 1:(nrow(IndMultiRepSize))) {
    for(j in 1:(nrow(IndMultiRepSize))) {
      Ind1 <- subset(datMultiRep.inc.perm.MultiSET, Subject==IndMultiRepSize$Subject[i])
      Ind2 <- subset(datMultiRep.inc.perm.MultiSET, Subject==IndMultiRepSize$Subject[j])
      if(length(table(Ind1$Setting))>1) { print("Ohoh, better check my data again!") }
      if(length(table(Ind2$Setting))>1) { print("Ohoh, better check my data again!") }
      if(Ind1$Setting[1]=="Captive" & Ind2$Setting[1]=="Captive") { SETwth.btw.Multi2 <- "Captive" }
      if(Ind1$Setting[1]=="Sanctuary" & Ind2$Setting[1]=="Sanctuary") { SETwth.btw.Multi2 <- "Sanctuary" }
      if(Ind1$Setting[1]=="Wild" & Ind2$Setting[1]=="Wild") { SETwth.btw.Multi2 <- "Wild" }
      SETdat.wth.btw.Multi2[which(rownames(SETdat.wth.btw.Multi2)==IndMultiRepSize$Subject[i]),which(colnames(SETdat.wth.btw.Multi2)==IndMultiRepSize$Subject[j])] <- SETwth.btw.Multi2
    }
  }
  ## set diagonal to NA (step6)
  SETwth.btw.Multi2 <- as.matrix(SETdat.wth.btw.Multi2)
  diag(SETwth.btw.Multi2) <- NA
  DiceMulti.C <- ifelse(SETwth.btw.Multi2=="Captive",DiceMulti,NA)
  DiceMulti.S <- ifelse(SETwth.btw.Multi2=="Sanctuary",DiceMulti,NA)
  DiceMulti.W <- ifelse(SETwth.btw.Multi2=="Wild",DiceMulti,NA)
  MDiceMultiC <- mean(DiceMulti.C[lower.tri(DiceMulti.C, diag = FALSE)],na.rm=TRUE)
  MDiceMultiS <- mean(DiceMulti.S[lower.tri(DiceMulti.S, diag = FALSE)],na.rm=TRUE)
  MDiceMultiW <- mean(DiceMulti.W[lower.tri(DiceMulti.W, diag = FALSE)],na.rm=TRUE)
  # omnibus statistic
  muMRS <- c(MDiceMultiC, MDiceMultiS, MDiceMultiW)
  mean.wMRS <- weighted.mean(muMRS, wSET=wSET, na.rm=TRUE)
  dat.Perm.MultiSET2[k] <- sum(wSET * (muMRS - mean.wMRS)^2, na.rm=TRUE)
  # pairwise differences
  perm.pairwiseMRS[k, ] <- c(
    MDiceMultiC - MDiceMultiS,
    MDiceMultiC - MDiceMultiW,
    MDiceMultiS - MDiceMultiW
  )
  flush.console()
  if(k %% 10 == 0) { print(paste0("Finished ", k, " out of ", nPerm, " simulations")) } 
}
hist(dat.Perm.MultiSET2)
abline(v=SS.empMRS, col="red") 

p.valueMRS <- (sum(dat.Perm.MultiSET2 >= SS.empMRS) + 1) / (length(dat.Perm.MultiSET2) + 1) #P-value = 0.000999001

#only if the variance differs:
pairwiseMRS.emp <- c(
  Emp.MDiceMultiC - Emp.MDiceMultiS,
  Emp.MDiceMultiC - Emp.MDiceMultiW,
  Emp.MDiceMultiS - Emp.MDiceMultiW
)

pvalsMRS <- sapply(1:3, function(i) {
  valid <- !is.na(perm.pairwiseMRS[, i])
  (sum(abs(perm.pairwiseMRS[valid, i]) >= abs(pairwiseMRS.emp[i])) + 1) /
    (sum(valid) + 1)
})

pvalsMRS.holm <- p.adjust(pvalsMRS, method = "holm")

contrast.namesSET <- c(
  "C-S","C-W","S-W")

resultsMRS <- data.frame(
  Contrast = contrast.namesSET,
  p_raw = pvalsMRS,
  p_holm = pvalsMRS.holm
)

resultsMRS
#  Contrast       p_raw      p_holm
#1      C-S 0.000999001 0.002997003 #
#2      C-W 0.000999001 0.002997003 #
#3      S-W 0.018981019 0.018981019 #


# Boxplot DiceMulti within/between (steps a)

library(ggplot2)
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
                         axis.title.y = element_text(size = 19, vjust = 2, family="sans"),
                         legend.text=  element_text(size = 15, family="sans", margin = margin(t = 10)),
                         legend.title = element_text(size = 15, vjust = 2, family="sans"),
                         legend.key = element_blank(),
                         legend.position = "right",
                         legend.spacing.x = unit(0.2, 'cm'),
                         title = element_text(size = 20, family="sans"),
                         strip.text = element_text(size = 15))

levels(WBddMultiSET$Setting.ind1) <- c("Captive", "Sanctuary", "Wild") # DiceMulti = dd ?

WBddMultiSETWithin <- subset(WBddMultiSET, WBddMultiSET$WithinBetween=="WithinS")

F2 <- ggplot() + 
  geom_boxplot(WBddMultiSETWithin, mapping=aes(x = Setting.ind1, y = DiceMulti), width = 0.9, fill = c("#068591", "#ba8211", "#ba0f09")) +
  geom_boxplot(WBddMultiSET, mapping=aes(x = WithinBetween, y = DiceMulti), width = 0.9, fill = c("grey90","grey90")) +
  geom_point(WBddMultiSETWithin, mapping=aes(x = Setting.ind1, y = DiceMulti), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_point(WBddMultiSET, mapping=aes(x = WithinBetween, y = DiceMulti), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey40") + 
  theme_angele_ss +
  ggtitle("(b)") +
  scale_y_continuous("Repertoire similarity among individuals") +
  scale_x_discrete(" ",
                   limits = c("Captive", "Sanctuary", "Wild", "BetweenS", "WithinS"),
                   labels = c("Within \nCaptive", "Within \nSanctuary", "Within \nWild", "Between \nsettings", "Within \nsettings"))+
  #facet_wrap(~DiceMulti, scales='free_x')+
  stat_summary(WBddMultiSETWithin, mapping=aes(x = Setting.ind1, y = DiceMulti), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3) +
  stat_summary(WBddMultiSET, mapping=aes(x = WithinBetween, y = DiceMulti), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

print(F2)

ggarrange(F1, F2, widths = c(1.7,1))

