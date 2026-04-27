4. GROUP MULTICOMPONENT REPERTOIRE SIMILARITY 

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


# Taux de similarit? entre les r?pertoires individuels (A-A/B-B, A-B): 

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

# step 2: define within- and between-group dyads 
## create data frame with Subject as rows and columns
dat.wth.btw.Multi <- data.frame(matrix(rep(NA, nrow(IndMultiRepSize) * nrow(IndMultiRepSize)), ncol=nrow(IndMultiRepSize)))
colnames(dat.wth.btw.Multi) <- IndMultiRepSize$Subject
rownames(dat.wth.btw.Multi) <- IndMultiRepSize$Subject
## loop over all pairs of individuals
for(i in 1:(nrow(IndMultiRepSize))) {
  for(j in 1:(nrow(IndMultiRepSize))) {
    Ind1 <- subset(datMultiRep, Subject==IndMultiRepSize$Subject[i])
    Ind2 <- subset(datMultiRep, Subject==IndMultiRepSize$Subject[j])
    if(length(table(Ind1$Group))>1) { print("Ohoh, better check my data again!") }
    if(length(table(Ind2$Group))>1) { print("Ohoh, better check my data again!") }
    if(Ind1$Group[1]==Ind2$Group[1]) { wth.btw.Multi <- "Within" } else { wth.btw.Multi <- "Between" }
    dat.wth.btw.Multi[which(rownames(dat.wth.btw.Multi)==IndMultiRepSize$Subject[i]),which(colnames(dat.wth.btw.Multi)==IndMultiRepSize$Subject[j])] <- wth.btw.Multi
  }
}
## set diagonal to NA
wth.btw.Multi <- as.matrix(dat.wth.btw.Multi)
diag(wth.btw.Multi) <- NA
#wth.btw.Multi <- as.data.frame(wth.btw.Multi)

# step 3: select within and between group dyads from the DiceMulti data frame (N=19)
DiceMulti.Within <- ifelse(wth.btw.Multi=="Within",DiceMulti,NA)
DiceMulti.Between <- ifelse(wth.btw.Multi=="Between",DiceMulti,NA)

(Emp.MDiceMultiWithin <- mean(DiceMulti.Within[lower.tri(DiceMulti.Within, diag = FALSE)],na.rm=TRUE)) # mean DiceMulti coefficient for within-dyads (= mean repertoire overlap between individuals in the same group (within, A-A and B-B))
# Mean DiceMulti within (A-A, B-B ...) groups (groupes A/B/C2/C4/L confondus) = 0.09810643
(Emp.MDiceMultiBetween <- mean(DiceMulti.Between[lower.tri(DiceMulti.Between, diag = FALSE)],na.rm=TRUE)) # mean DiceMulti coefficient for between-dyads (= mean repertoire overlap between individuals in different groups (between, A-B))
# Mean DiceMulti between (A-B, A-C2, A-C4 ...) groups = 0.04618709


# Repertoire similarity between groups Within et Between

# step 4: Combine all info in 1 table (1 line = 1 dyad --> create data frame dd with columns Ind1, Ind2, Dice, Within/Between)
## Tranform Within matrix in table
yWithinMulti <- expand.grid(rownames(DiceMulti.Within), colnames(DiceMulti.Within)) 
labsWithin <- yWithinMulti[as.vector(upper.tri(DiceMulti.Within, diag = FALSE)), ] 
yWithinMulti <- cbind(labsWithin, DiceMulti.Within[upper.tri(DiceMulti.Within,diag=FALSE)]) 
colnames(yWithinMulti) <- c("Ind1","Ind2","DiceMulti")
yWithinMulti <- yWithinMulti[!is.na(yWithinMulti$DiceMulti), ] 
yWithinMulti$WithinBetween <- rep("Within",nrow(yWithinMulti)) 
mean(yWithinMulti$Dice)

## Tranform Between matrix in table
yBetweenMulti <- expand.grid(rownames(DiceMulti.Between), colnames(DiceMulti.Between))
labsBetween <- yBetweenMulti[as.vector(upper.tri(DiceMulti.Between, diag = FALSE)), ]
yBetweenMulti <- cbind(labsBetween, DiceMulti.Between[upper.tri(DiceMulti.Between,diag=FALSE)])
colnames(yBetweenMulti) <- c("Ind1","Ind2","DiceMulti")
yBetweenMulti <- yBetweenMulti[!is.na(yBetweenMulti$DiceMulti), ]
yBetweenMulti$WithinBetween <- rep("Between",nrow(yBetweenMulti))
mean(yBetweenMulti$Dice)

##Combine Within and Between tables
WBddMulti <- rbind(yWithinMulti,yBetweenMulti)

datMultiRep.sub <- datMultiRep[,c("Subject","Group")] 
datMultiRep.sub <- datMultiRep.sub[!duplicated(datMultiRep.sub), ] 
WBddMulti <- merge(WBddMulti, datMultiRep.sub, by.x="Ind1", by.y="Subject", sort=FALSE, all.y=FALSE) 
WBddMulti <- merge(WBddMulti, datMultiRep.sub, by.x="Ind2", by.y="Subject", sort=FALSE) 
colnames(WBddMulti)[which(colnames(WBddMulti)=="Group.x")] <- "Group.ind1"
colnames(WBddMulti)[which(colnames(WBddMulti)=="Group.y")] <- "Group.ind2"
# V?rifier l'ordre des colonnes
WBddMulti<-WBddMulti[,c("Ind1","Ind2","DiceMulti","WithinBetween","Group.ind1","Group.ind2")]

# step 5: matrix permutation test ---
nPerm <- 1000
dat.Perm.Multi <- c()
for(k in 1:nPerm) {
  ## permute group in datMultiRep
  ind.group.perm.Multi <- datMultiRep.sub
  ind.group.perm.Multi$Group <- sample(ind.group.perm.Multi$Group, nrow(ind.group.perm.Multi), replace=FALSE)
  datMultiRep.inc.perm.Multi <- datMultiRep[,-which(colnames(datMultiRep)=="Group")] 
  datMultiRep.inc.perm.Multi <- merge(datMultiRep.inc.perm.Multi, ind.group.perm.Multi, by.x="Subject", by.y="Subject", sort=FALSE)
  ## repeat step 2 & 3
  ## create data frame with Subject as rows and columns (step2)
  dat.wth.btw.Multi <- data.frame(matrix(rep(NA, nrow(IndMultiRepSize) * nrow(IndMultiRepSize)), ncol=nrow(IndMultiRepSize)))
  colnames(dat.wth.btw.Multi) <- IndMultiRepSize$Subject
  rownames(dat.wth.btw.Multi) <- IndMultiRepSize$Subject
  ## loop over all pairs of individuals (step2)
  for(i in 1:(nrow(IndMultiRepSize))) {
    for(j in 1:(nrow(IndMultiRepSize))) {
      Ind1 <- subset(datMultiRep.inc.perm.Multi, Subject==IndMultiRepSize$Subject[i])
      Ind2 <- subset(datMultiRep.inc.perm.Multi, Subject==IndMultiRepSize$Subject[j])
      if(length(table(Ind1$Group))>1) { print("Ohoh, better check my data again!") }
      if(length(table(Ind2$Group))>1) { print("Ohoh, better check my data again!") }
      if(Ind1$Group[1]==Ind2$Group[1]) { wth.btw.Multi <- "Within" } else { wth.btw.Multi <- "Between" }
      dat.wth.btw.Multi[which(rownames(dat.wth.btw.Multi)==IndMultiRepSize$Subject[i]),which(colnames(dat.wth.btw.Multi)==IndMultiRepSize$Subject[j])] <- wth.btw.Multi
    }
  }
  ## set diagonal to NA (step3)
  wth.btw.Multi <- as.matrix(dat.wth.btw.Multi)
  diag(wth.btw.Multi) <- NA
  DiceMulti.Within <- ifelse(wth.btw.Multi=="Within",DiceMulti,NA)
  DiceMulti.Between <- ifelse(wth.btw.Multi=="Between",DiceMulti,NA)
  MDiceMultiWithin <- mean(DiceMulti.Within[lower.tri(DiceMulti.Within, diag = FALSE)],na.rm=TRUE)
  MDiceMultiBetween <- mean(DiceMulti.Between[lower.tri(DiceMulti.Between, diag = FALSE)],na.rm=TRUE)
  dat.Perm.Multi[k] <- MDiceMultiWithin - MDiceMultiBetween 
  flush.console()
  if(k %% 10 == 0) { print(paste0("Finished ", k, " out of ", nPerm, " simulations")) } 
}

hist(dat.Perm.Multi)
abline(v=Emp.MDiceMultiWithin - Emp.MDiceMultiBetween, col="red") 
pMRG <- (sum(abs(dat.Perm.Multi) >= abs(Emp.MDiceMultiWithin - Emp.MDiceMultiBetween)) + 1) / (nPerm + 1)
# P-value = 0.000999


# Significance thresholds are P>=0.975 and P<=0.025. 
# This is because we are looking at the deviation from 0, which can either be negative or positive.
# Because the distribution of differences is not necessarily symmetric around zero, we cannot calculate P-values using absolute values.


# Repertoire similarity within group A (A-A), within group B (B-B), within group C2, within group C4, within group L, within group K

# step 2b: define within and between setting diads --- for A, B, C2, C4, L and K groups 
## create data frame with Subject as rows and columns
dat.wth.btw.Multi2 <- data.frame(matrix(rep(NA, nrow(IndMultiRepSize) * nrow(IndMultiRepSize)), ncol=nrow(IndMultiRepSize)))
colnames(dat.wth.btw.Multi2) <- IndMultiRepSize$Subject
rownames(dat.wth.btw.Multi2) <- IndMultiRepSize$Subject
## loop over all pairs of individuals
for(i in 1:(nrow(IndMultiRepSize))) {
  for(j in 1:(nrow(IndMultiRepSize))) {
    wth.btw.Multi2 <- NA
    Ind1 <- subset(datMultiRep, Subject==IndMultiRepSize$Subject[i])
    Ind2 <- subset(datMultiRep, Subject==IndMultiRepSize$Subject[j])
    if(length(table(Ind1$Group))>1) { print("Ohoh, better check my data again!") }
    if(length(table(Ind2$Group))>1) { print("Ohoh, better check my data again!") }
    if(Ind1$Group[1]=="A" & Ind2$Group[1]=="A") { wth.btw.Multi2 <- "A" }
    if(Ind1$Group[1]=="B" & Ind2$Group[1]=="B") { wth.btw.Multi2 <- "B" }
    if(Ind1$Group[1]=="C2" & Ind2$Group[1]=="C2") { wth.btw.Multi2 <- "C2" }
    if(Ind1$Group[1]=="C4" & Ind2$Group[1]=="C4") { wth.btw.Multi2 <- "C4" }
    if(Ind1$Group[1]=="L" & Ind2$Group[1]=="L") { wth.btw.Multi2 <- "L" }
    if(Ind1$Group[1]=="K" & Ind2$Group[1]=="K") { wth.btw.Multi2 <- "K" }
    dat.wth.btw.Multi2[which(rownames(dat.wth.btw.Multi2)==IndMultiRepSize$Subject[i]),which(colnames(dat.wth.btw.Multi2)==IndMultiRepSize$Subject[j])] <- wth.btw.Multi2
  }
}
## set diagonal to NA
wth.btw.Multi2 <- as.matrix(dat.wth.btw.Multi2)
diag(wth.btw.Multi2) <- NA
#wth.btw.Multi2 <- as.data.frame(wth.btw.Multi2)

# step 3b: select within and between group diads from the DiceMulti data frame ---
DiceMulti.A <- ifelse(wth.btw.Multi2=="A",DiceMulti,NA)
DiceMulti.B <- ifelse(wth.btw.Multi2=="B",DiceMulti,NA)
DiceMulti.C2 <- ifelse(wth.btw.Multi2=="C2",DiceMulti,NA)
DiceMulti.C4 <- ifelse(wth.btw.Multi2=="C4",DiceMulti,NA)
DiceMulti.L <- ifelse(wth.btw.Multi2=="L",DiceMulti,NA)
DiceMulti.K <- ifelse(wth.btw.Multi2=="K",DiceMulti,NA)

(Emp.MDiceMultiA <- mean(DiceMulti.A[lower.tri(DiceMulti.A, diag = FALSE)],na.rm=TRUE))
# Mean DiceMulti within group A = 0.1389911 --> big individual variation (0 =  no overlap ; 1 = 100% overlap)
(Emp.MDiceMultiB <- mean(DiceMulti.B[lower.tri(DiceMulti.B, diag = FALSE)],na.rm=TRUE))
# Mean DiceMulti within group B = 0.2857143
(Emp.MDiceMultiC2 <- mean(DiceMulti.C2[lower.tri(DiceMulti.C2, diag = FALSE)],na.rm=TRUE))
# Mean DiceMulti within group C2 = 0.1099815
(Emp.MDiceMultiC4 <- mean(DiceMulti.C4[lower.tri(DiceMulti.C4, diag = FALSE)],na.rm=TRUE))
# Mean DiceMulti within group C4 = 0
(Emp.MDiceMultiL <- mean(DiceMulti.L[lower.tri(DiceMulti.L, diag = FALSE)],na.rm=TRUE))
# Mean DiceMulti within group L = 0.05314279
(Emp.MDiceMultiK <- mean(DiceMulti.K[lower.tri(DiceMulti.K, diag = FALSE)],na.rm=TRUE))
# Mean DiceMulti within group K = 0.1055461

## step 4b: Combiner toutes les infos dans un tableau (1 ligne = 1 diade --> create data frame ddMulti with columns Ind1, Ind2, DiceMulti, A/B)
## Tranform A matrix in table
yAMulti <- expand.grid(rownames(DiceMulti.A), colnames(DiceMulti.A)) 
labsA <- yAMulti[as.vector(upper.tri(DiceMulti.A, diag = FALSE)), ] 
yAMulti <- cbind(labsA, DiceMulti.A[upper.tri(DiceMulti.A,diag=FALSE)]) 
colnames(yAMulti) <- c("Ind1","Ind2","DiceMulti")
yAMulti <- yAMulti[!is.na(yAMulti$DiceMulti), ] 
yAMulti$Group <- rep("A",nrow(yAMulti)) 
mean(yAMulti$Dice)

## Tranform B matrix in table
yBMulti <- expand.grid(rownames(DiceMulti.B), colnames(DiceMulti.B))
labsB <- yBMulti[as.vector(upper.tri(DiceMulti.B, diag = FALSE)), ]
yBMulti <- cbind(labsB, DiceMulti.B[upper.tri(DiceMulti.B,diag=FALSE)])
colnames(yBMulti) <- c("Ind1","Ind2","DiceMulti")
yBMulti <- yBMulti[!is.na(yBMulti$DiceMulti), ]
yBMulti$Group <- rep("B",nrow(yBMulti))
mean(yBMulti$Dice)

## Tranform C2 matrix in table
yC2Multi <- expand.grid(rownames(DiceMulti.C2), colnames(DiceMulti.C2)) 
labsC2 <- yC2Multi[as.vector(upper.tri(DiceMulti.C2, diag = FALSE)), ] 
yC2Multi <- cbind(labsC2, DiceMulti.C2[upper.tri(DiceMulti.C2,diag=FALSE)]) 
colnames(yC2Multi) <- c("Ind1","Ind2","DiceMulti")
yC2Multi <- yC2Multi[!is.na(yC2Multi$DiceMulti), ] 
yC2Multi$Group <- rep("C2",nrow(yC2Multi)) 
mean(yC2Multi$Dice)

## Tranform C4 matrix in table
yC4Multi <- expand.grid(rownames(DiceMulti.C4), colnames(DiceMulti.C4)) 
labsC4 <- yC4Multi[as.vector(upper.tri(DiceMulti.C4, diag = FALSE)), ] 
yC4Multi <- cbind(labsC4, DiceMulti.C4[upper.tri(DiceMulti.C4,diag=FALSE)]) 
colnames(yC4Multi) <- c("Ind1","Ind2","DiceMulti")
yC4Multi <- yC4Multi[!is.na(yC4Multi$DiceMulti), ] 
yC4Multi$Group <- rep("C4",nrow(yC4Multi)) 
mean(yC4Multi$Dice)

## Tranform L matrix in table
yLMulti <- expand.grid(rownames(DiceMulti.L), colnames(DiceMulti.L)) 
labsL <- yLMulti[as.vector(upper.tri(DiceMulti.L, diag = FALSE)), ] 
yLMulti <- cbind(labsL, DiceMulti.L[upper.tri(DiceMulti.L,diag=FALSE)]) 
colnames(yLMulti) <- c("Ind1","Ind2","DiceMulti")
yLMulti <- yLMulti[!is.na(yLMulti$DiceMulti), ] 
yLMulti$Group <- rep("L",nrow(yLMulti)) 
mean(yLMulti$Dice)

## Tranform K matrix in table
yKMulti <- expand.grid(rownames(DiceMulti.K), colnames(DiceMulti.K)) 
labsK <- yKMulti[as.vector(upper.tri(DiceMulti.K, diag = FALSE)), ] 
yKMulti <- cbind(labsK, DiceMulti.K[upper.tri(DiceMulti.K,diag=FALSE)]) 
colnames(yKMulti) <- c("Ind1","Ind2","DiceMulti")
yKMulti <- yKMulti[!is.na(yKMulti$DiceMulti), ] 
yKMulti$Group <- rep("K",nrow(yKMulti)) 
mean(yKMulti$Dice)

## Combine A/B/C2/C4/L/K tables
ddMulti <- rbind(yAMulti,yBMulti,yC2Multi,yC4Multi,yLMulti,yKMulti)

datMultiRep.sub <- datMultiRep[,c("Subject","Group")] 
datMultiRep.sub <- datMultiRep.sub[!duplicated(datMultiRep.sub), ] 
ddMulti <- merge(ddMulti, datMultiRep.sub, by.x="Ind1", by.y="Subject", sort=FALSE, all.y=FALSE) 
ddMulti <- merge(ddMulti, datMultiRep.sub, by.x="Ind2", by.y="Subject", sort=FALSE) 
colnames(ddMulti)[which(colnames(ddMulti)=="Group.x")] <- "Group.ind1"
colnames(ddMulti)[which(colnames(ddMulti)=="Group.y")] <- "Group.ind2"
# Order columns
ddMulti<-ddMulti[,c("Ind1","Ind2","DiceMulti","Group","Group.ind1","Group.ind2")]

mu.empMRG <- c(Emp.MDiceMultiA, Emp.MDiceMultiB, Emp.MDiceMultiC2, Emp.MDiceMultiC4, Emp.MDiceMultiL, Emp.MDiceMultiK)
names(mu.empMRG) <- c("A","B","C2","C4","L","K")
# Get the number of individuals per group
group_sizesGRPE <- table(datMultiRep.sub$Group)[c("A","B","C2","C4","L","K")]
wGRPE <- group_sizesGRPE * (group_sizesGRPE - 1) / 2
# Compute the weighted sums of squares
mean.w.empMRG <- weighted.mean(mu.empMRG, wGRPE=wGRPE, na.rm=TRUE)
SS.empMRG <- sum(wGRPE * (mu.empMRG - mean.w.empMRG)^2, na.rm=TRUE) 


# step 5b: matrix permutation test ---
nPerm <- 1000
dat.Perm.Multi2 <- c()
perm.pairwiseMRG <- matrix(NA_real_, nrow = nPerm, ncol = 15)
for(k in 1:nPerm) {
  ## permute group in datMultiRep
  ind.group.perm.Multi <- datMultiRep.sub
  ind.group.perm.Multi$Group <- sample(ind.group.perm.Multi$Group, nrow(ind.group.perm.Multi), replace=FALSE)
  datMultiRep.inc.perm.Multi <- datMultiRep[,-which(colnames(datMultiRep)=="Group")] 
  datMultiRep.inc.perm.Multi <- merge(datMultiRep.inc.perm.Multi, ind.group.perm.Multi, by.x="Subject", by.y="Subject", sort=FALSE)
  ## repeat step 2 & 3
  ## create data frame with Subject as rows and columns (step2)
  dat.wth.btw.Multi2 <- data.frame(matrix(rep(NA, nrow(IndMultiRepSize) * nrow(IndMultiRepSize)), ncol=nrow(IndMultiRepSize)))
  colnames(dat.wth.btw.Multi2) <- IndMultiRepSize$Subject
  rownames(dat.wth.btw.Multi2) <- IndMultiRepSize$Subject
  ## loop over all pairs of individuals (step2)
  for(i in 1:(nrow(IndMultiRepSize))) {
    for(j in 1:(nrow(IndMultiRepSize))) {
      Ind1 <- subset(datMultiRep.inc.perm.Multi, Subject==IndMultiRepSize$Subject[i])
      Ind2 <- subset(datMultiRep.inc.perm.Multi, Subject==IndMultiRepSize$Subject[j])
      if(length(table(Ind1$Group))>1) { print("Ohoh, better check my data again!") }
      if(length(table(Ind2$Group))>1) { print("Ohoh, better check my data again!") }
      if(Ind1$Group[1]=="A" & Ind2$Group[1]=="A") { wth.btw.Multi2 <- "A" }
      if(Ind1$Group[1]=="B" & Ind2$Group[1]=="B") { wth.btw.Multi2 <- "B" }
      if(Ind1$Group[1]=="C2" & Ind2$Group[1]=="C2") { wth.btw.Multi2 <- "C2" }
      if(Ind1$Group[1]=="C4" & Ind2$Group[1]=="C4") { wth.btw.Multi2 <- "C4" }
      if(Ind1$Group[1]=="L" & Ind2$Group[1]=="L") { wth.btw.Multi2 <- "L" }
      if(Ind1$Group[1]=="K" & Ind2$Group[1]=="K") { wth.btw.Multi2 <- "K" }
      dat.wth.btw.Multi2[which(rownames(dat.wth.btw.Multi2)==IndMultiRepSize$Subject[i]),which(colnames(dat.wth.btw.Multi2)==IndMultiRepSize$Subject[j])] <- wth.btw.Multi2
    }
  }
  ## set diagonal to NA (step6)
  wth.btw.Multi2 <- as.matrix(dat.wth.btw.Multi2)
  diag(wth.btw.Multi2) <- NA
  DiceMulti.A <- ifelse(wth.btw.Multi2=="A",DiceMulti,NA)
  DiceMulti.B <- ifelse(wth.btw.Multi2=="B",DiceMulti,NA)
  DiceMulti.C2 <- ifelse(wth.btw.Multi2=="C2",DiceMulti,NA)
  DiceMulti.C4 <- ifelse(wth.btw.Multi2=="C4",DiceMulti,NA)
  DiceMulti.L <- ifelse(wth.btw.Multi2=="L",DiceMulti,NA)
  DiceMulti.K <- ifelse(wth.btw.Multi2=="K",DiceMulti,NA)
  MDiceMultiA <- mean(DiceMulti.A[lower.tri(DiceMulti.A, diag = FALSE)],na.rm=TRUE)
  MDiceMultiB <- mean(DiceMulti.B[lower.tri(DiceMulti.B, diag = FALSE)],na.rm=TRUE)
  MDiceMultiC2 <- mean(DiceMulti.C2[lower.tri(DiceMulti.C2, diag = FALSE)],na.rm=TRUE)
  MDiceMultiC4 <- mean(DiceMulti.C4[lower.tri(DiceMulti.C4, diag = FALSE)],na.rm=TRUE)
  MDiceMultiL <- mean(DiceMulti.L[lower.tri(DiceMulti.L, diag = FALSE)],na.rm=TRUE)
  MDiceMultiK <- mean(DiceMulti.K[lower.tri(DiceMulti.K, diag = FALSE)],na.rm=TRUE)
  # omnibus statistic
  muMRG <- c(MDiceMultiA, MDiceMultiB, MDiceMultiC2, MDiceMultiC4, MDiceMultiL, MDiceMultiK)
  mean.wMRG <- weighted.mean(muMRG, wGRPE=wGRPE, na.rm=TRUE)
  dat.Perm.Multi2[k] <- sum(wGRPE * (muMRG - mean.wMRG)^2, na.rm=TRUE)
  # pairwise differences
  perm.pairwiseMRG[k, ] <- c(
    MDiceMultiA - MDiceMultiB,
    MDiceMultiA - MDiceMultiC2,
    MDiceMultiA - MDiceMultiC4,
    MDiceMultiA - MDiceMultiL,
    MDiceMultiA - MDiceMultiK,
    MDiceMultiB - MDiceMultiC2,
    MDiceMultiB - MDiceMultiC4,
    MDiceMultiB - MDiceMultiL,
    MDiceMultiB - MDiceMultiK,
    MDiceMultiC2 - MDiceMultiC4,
    MDiceMultiC2 - MDiceMultiL,
    MDiceMultiC2 - MDiceMultiK,
    MDiceMultiC4 - MDiceMultiL,
    MDiceMultiC4 - MDiceMultiK,
    MDiceMultiL - MDiceMultiK
  )
  flush.console()
  if(k %% 10 == 0) { print(paste0("Finished ", k, " out of ", nPerm, " simulations")) } 
}
hist(dat.Perm.Multi2)
abline(v=SS.empMRG, col="red") 

p.valueMRG <- (sum(dat.Perm.Multi2 >= SS.empMRG) + 1) / (length(dat.Perm.Multi2) + 1) #P-value = 0.000999

#only if the variance differs:
pairwiseMRG.emp <- c(
  Emp.MDiceMultiA - Emp.MDiceMultiB,
  Emp.MDiceMultiA - Emp.MDiceMultiC2,
  Emp.MDiceMultiA - Emp.MDiceMultiC4,
  Emp.MDiceMultiA - Emp.MDiceMultiL,
  Emp.MDiceMultiA - Emp.MDiceMultiK,
  Emp.MDiceMultiB - Emp.MDiceMultiC2,
  Emp.MDiceMultiB - Emp.MDiceMultiC4,
  Emp.MDiceMultiB - Emp.MDiceMultiL,
  Emp.MDiceMultiB - Emp.MDiceMultiK,
  Emp.MDiceMultiC2 - Emp.MDiceMultiC4,
  Emp.MDiceMultiC2 - Emp.MDiceMultiL,
  Emp.MDiceMultiC2 - Emp.MDiceMultiK,
  Emp.MDiceMultiC4 - Emp.MDiceMultiL,
  Emp.MDiceMultiC4 - Emp.MDiceMultiK,
  Emp.MDiceMultiL - Emp.MDiceMultiK
)

pvalsMRG <- sapply(1:15, function(i) {
  valid <- !is.na(perm.pairwiseMRG[, i])
  (sum(abs(perm.pairwiseMRG[valid, i]) >= abs(pairwiseMRG.emp[i])) + 1) /
    (sum(valid) + 1)
})

pvalsMRG.holm <- p.adjust(pvalsMRG, method = "holm")

contrast.namesGRPE <- c(
  "A-B","A-C2","A-C4","A-L","A-K",
  "B-C2","B-C4","B-L","B-K",
  "C2-C4","C2-L","C2-K",
  "C4-L","C4-K",
  "L-K"
)

resultsMRG <- data.frame(
  Contrast = contrast.namesGRPE,
  p_raw = pvalsMRG,
  p_holm = pvalsMRG.holm
)

resultsMRG
#   Contrast       p_raw     p_holm
#1       A-B 0.000999001 0.01498501 #
#2      A-C2 0.005994006 0.01498501 #
#3      A-C4 0.000999001 0.01498501 #
#4       A-L 0.000999001 0.01498501 #
#5       A-K 0.002997003 0.01498501 #
#6      B-C2 0.000999001 0.01498501 #
#7      B-C4 0.000999001 0.01498501 #
#8       B-L 0.000999001 0.01498501 #
#9       B-K 0.000999001 0.01498501 #
#10    C2-C4 0.000999001 0.01498501 #
#11     C2-L 0.000999001 0.01498501 #
#12     C2-K 0.610389610 0.61038961
#13     C4-L 0.000999001 0.01498501 #
#14     C4-K 0.000999001 0.01498501 #
#15      L-K 0.000999001 0.01498501 #


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
                         axis.title.y = element_text(size = 15, vjust = 2, family="sans"),
                         legend.text=  element_text(size = 15, family="sans", margin = margin(t = 10)),
                         legend.title = element_text(size = 15, vjust = 2, family="sans"),
                         legend.key = element_blank(),
                         legend.position = "right",
                         legend.spacing.x = unit(0.2, 'cm'),
                         title = element_text(size = 25, family="sans"),
                         strip.text = element_text(size = 15))

levels(WBddMulti$Group.ind1) <- c("A", "B", "L", "C2", "C4", "K") # dice = ddMulti ?

levels(ddMulti$Group.ind1) <- c("A", "B", "L", "C2", "C4", "K") # dice = ddMulti ?

WBddMultiWithin <- subset(WBddMulti, WBddMulti$WithinBetween=="Within")

F1 <- ggplot() + 
  geom_boxplot(WBddMultiWithin, mapping=aes(x = Group.ind1, y = DiceMulti), width = 0.9, fill = c("#007ABB", "#00AFBB","#259C39", "#E7B800","#E79A00", "#ba0f09")) +
  geom_boxplot(WBddMulti, mapping=aes(x = WithinBetween, y = DiceMulti), width = 0.9, fill = c("grey90","grey90")) +
  geom_point(WBddMultiWithin, mapping=aes(x = Group.ind1, y = DiceMulti), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_point(WBddMulti, mapping=aes(x = WithinBetween, y = DiceMulti), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "grey40") + 
  theme_angele_ss +
  ggtitle("(a)") +
  scale_y_continuous("Repertoire similarity among individuals") +
  scale_x_discrete(" ",
                   limits = c("A", "B", "L", "C2", "C4", "K", "Between", "Within"),
                   labels = c("Within \ngroup A", "Within \ngroup B", "Within \ngroup L", "Within \ngroup C2", "Within \ngroup C4", "Within \ngroup K", "Between \ngroups", "Within \ngroups"))+
  #facet_wrap(~DiceMulti, scales='free_x')+
  stat_summary(WBddMultiWithin, mapping=aes(x = Group.ind1, y = DiceMulti), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3) +
  stat_summary(WBddMulti, mapping=aes(x = WithinBetween, y = DiceMulti), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3)

# + theme(axis.text.x = element_blank()) 

print(F1)
