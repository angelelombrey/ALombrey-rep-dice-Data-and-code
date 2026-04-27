4. GROUP SIGNAL REPERTOIRE SIMILARITY 

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
# Create a "Dyad" column ensuring that dyads are treated identically
Data$Dyad <- apply(Data[, c("Subject", "Recipient")], 1, function(row) {
  paste(sort(row), collapse = "/")
})

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

# Omit all levels of Subject that contributed fewer than 30 cases --- 118 --> 93 individuals
nb.signals.per.ind <- data.frame(Data %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(nb.signals.per.ind) <- c("Subject","NsignalsTot") 
nb.signals.per.ind$Subject <- as.character(nb.signals.per.ind$Subject)
nb.signals.per.ind <- subset(nb.signals.per.ind, nb.signals.per.ind$NsignalsTot>=30) # Keeps only the values in "column" that appear more than 30 times (N >= 30)
Data <- subset(Data, Subject %in% nb.signals.per.ind$Subject)

Data <- Data[!is.na(Data$Group), ]

IndRepSize <- aggregate(data=Data, Behavior ~ Subject, function(x) length(unique(x))) ; colnames(IndRepSize) <- c("Subject","RepertoireSize") 
dat.ind.bhv2 <- data.frame(Data %>% group_by(Subject,Behavior) %>% dplyr::summarize(count=n())) ; colnames(dat.ind.bhv2) <- c("Subject","Behavior","NsignalsPbehavior") 

________________________________________________________________________________________________________________________________________________________________________________________________


# Taux de similarit? entre les r?pertoires individuels (A-A/B-B, A-B): 

# Step 1 : calculate Dice coefficient based on formula: dc = (2 x number of behaviours two inds have in common)/(R of ind1 + R ind2) ---
## create data frame in which every combination of individual and behaviour is counted
dat.ind.bhv <- data.frame(table(Data$Subject,Data$Behavior)) ; colnames(dat.ind.bhv) <- c("Subject","Behavior","Nsignals")
## create data frame with Subject as rows and columns
dat.ind.ind <- data.frame(matrix(rep(NA, nrow(IndRepSize) * nrow(IndRepSize)), ncol=nrow(IndRepSize)))
colnames(dat.ind.ind) <- IndRepSize$Subject
rownames(dat.ind.ind) <- IndRepSize$Subject
## loop over all pairs of individuals
for(i in 1:(nrow(IndRepSize))) {
  for(j in 1:(nrow(IndRepSize))) {
    Ind1 <- subset(dat.ind.bhv2, Subject==IndRepSize$Subject[i])
    Ind2 <- subset(dat.ind.bhv2, Subject==IndRepSize$Subject[j])
    ovrlp <- length(which(Ind1$Behavior %in% Ind2$Behavior)) # number of behaviors 2 individuals have in common
    RInd1 <- subset(IndRepSize, Subject==IndRepSize$Subject[i])$RepertoireSize	# which is the same as nrow(Ind1)
    RInd2 <- subset(IndRepSize, Subject==IndRepSize$Subject[j])$RepertoireSize	# which is the same as nrow(Ind2)
    dat.ind.ind[which(rownames(dat.ind.ind)==IndRepSize$Subject[i]),which(colnames(dat.ind.ind)==IndRepSize$Subject[j])] <- (2 * ovrlp) / (RInd1 + RInd2)
  }
}
## set diagonal to NA
Dice <- as.matrix(dat.ind.ind)
diag(Dice) <- NA
#Dice <- as.data.frame(Dice)

# step 2: define within- and between-group dyads 
## create data frame with Subject as rows and columns
dat.wth.btw <- data.frame(matrix(rep(NA, nrow(IndRepSize) * nrow(IndRepSize)), ncol=nrow(IndRepSize)))
colnames(dat.wth.btw) <- IndRepSize$Subject
rownames(dat.wth.btw) <- IndRepSize$Subject
## loop over all pairs of individuals
for(i in 1:(nrow(IndRepSize))) {
  for(j in 1:(nrow(IndRepSize))) {
    Ind1 <- subset(Data, Subject==IndRepSize$Subject[i])
    Ind2 <- subset(Data, Subject==IndRepSize$Subject[j])
    if(length(table(Ind1$Group))>1) { print("Ohoh, better check my data again!") }
    if(length(table(Ind2$Group))>1) { print("Ohoh, better check my data again!") }
    if(Ind1$Group[1]==Ind2$Group[1]) { wth.btw <- "Within" } else { wth.btw <- "Between" }
    dat.wth.btw[which(rownames(dat.wth.btw)==IndRepSize$Subject[i]),which(colnames(dat.wth.btw)==IndRepSize$Subject[j])] <- wth.btw
  }
}
## set diagonal to NA
wth.btw <- as.matrix(dat.wth.btw)
diag(wth.btw) <- NA
#wth.btw <- as.data.frame(wth.btw)

# step 3: select within and between group dyads from the Dice data frame (N=19)
Dice.Within <- ifelse(wth.btw=="Within",Dice,NA)
Dice.Between <- ifelse(wth.btw=="Between",Dice,NA)

(Emp.MDiceWithin <- mean(Dice.Within[lower.tri(Dice.Within, diag = FALSE)],na.rm=TRUE)) # mean Dice coefficient for within-dyads (= mean repertoire overlap between individuals in the same group (within, A-A and B-B))
# Mean Dice within (A-A, B-B ...) groups (groupes A/B/C2/C4/L confondus) = 0.5543321
(Emp.MDiceBetween <- mean(Dice.Between[lower.tri(Dice.Between, diag = FALSE)],na.rm=TRUE)) # mean Dice coefficient for between-dyads (= mean repertoire overlap between individuals in different groups (between, A-B))
# Mean Dice between (A-B, A-C2, A-C4 ...) groups = 0.4679664

# Repertoire similarity entre les groupes Within et Between

# step 4: Combine all info in 1 table (1 line = 1 dyad --> create data frame dd with columns Ind1, Ind2, Dice, Within/Between)
## Tranform Within matrix in table
yWithin <- expand.grid(rownames(Dice.Within), colnames(Dice.Within)) # associate each line name (= Ind1) to the corresponding column name (=Ind2) 
labsWithin <- yWithin[as.vector(upper.tri(Dice.Within, diag = FALSE)), ] # take out lower half of the matrix (=duplicates) and the diagonal
yWithin <- cbind(labsWithin, Dice.Within[upper.tri(Dice.Within,diag=FALSE)]) # add one column with the Dice coef for each dyad
colnames(yWithin) <- c("Ind1","Ind2","Dice")
yWithin <- yWithin[!is.na(yWithin$Dice), ] # take out Dice = NA
yWithin$WithinBetween <- rep("Within",nrow(yWithin)) # add column Within

## Tranform Between matrix in table
yBetween <- expand.grid(rownames(Dice.Between), colnames(Dice.Between))
labsBetween <- yBetween[as.vector(upper.tri(Dice.Between, diag = FALSE)), ]
yBetween <- cbind(labsBetween, Dice.Between[upper.tri(Dice.Between,diag=FALSE)])
colnames(yBetween) <- c("Ind1","Ind2","Dice")
yBetween <- yBetween[!is.na(yBetween$Dice), ]
yBetween$WithinBetween <- rep("Between",nrow(yBetween))

##Combine Within and Between tables
WBdd <- rbind(yWithin,yBetween)

Data.sub <- Data[,c("Subject","Group")] 
Data.sub <- Data.sub[!duplicated(Data.sub), ] 
WBdd <- merge(WBdd, Data.sub, by.x="Ind1", by.y="Subject", sort=FALSE, all.y=FALSE) 
WBdd <- merge(WBdd, Data.sub, by.x="Ind2", by.y="Subject", sort=FALSE) 
colnames(WBdd)[which(colnames(WBdd)=="Group.x")] <- "Group.ind1"
colnames(WBdd)[which(colnames(WBdd)=="Group.y")] <- "Group.ind2"
# Order columns
WBdd<-WBdd[,c("Ind1","Ind2","Dice","WithinBetween","Group.ind1","Group.ind2")]

# step 5: matrix permutation test ---
nPerm <- 1000
dat.Perm <- c()
for(k in 1:nPerm) {
  ## permute group in Data
  ind.group.perm <- Data.sub
  ind.group.perm$Group <- sample(ind.group.perm$Group, nrow(ind.group.perm), replace=FALSE)
  Data.inc.perm <- Data[,-which(colnames(Data)=="Group")] 
  Data.inc.perm <- merge(Data.inc.perm, ind.group.perm, by.x="Subject", by.y="Subject", sort=FALSE)
  ## repeat step 2 & 3
  ## create data frame with Subject as rows and columns (step2)
  dat.wth.btw <- data.frame(matrix(rep(NA, nrow(IndRepSize) * nrow(IndRepSize)), ncol=nrow(IndRepSize)))
  colnames(dat.wth.btw) <- IndRepSize$Subject
  rownames(dat.wth.btw) <- IndRepSize$Subject
  ## loop over all pairs of individuals (step2)
  for(i in 1:(nrow(IndRepSize))) {
    for(j in 1:(nrow(IndRepSize))) {
      Ind1 <- subset(Data.inc.perm, Subject==IndRepSize$Subject[i])
      Ind2 <- subset(Data.inc.perm, Subject==IndRepSize$Subject[j])
      if(length(table(Ind1$Group))>1) { print("Ohoh, better check my data again!") }
      if(length(table(Ind2$Group))>1) { print("Ohoh, better check my data again!") }
      if(Ind1$Group[1]==Ind2$Group[1]) { wth.btw <- "Within" } else { wth.btw <- "Between" }
      dat.wth.btw[which(rownames(dat.wth.btw)==IndRepSize$Subject[i]),which(colnames(dat.wth.btw)==IndRepSize$Subject[j])] <- wth.btw
    }
  }
  ## set diagonal to NA (step3)
  wth.btw <- as.matrix(dat.wth.btw)
  diag(wth.btw) <- NA
  Dice.Within <- ifelse(wth.btw=="Within",Dice,NA)
  Dice.Between <- ifelse(wth.btw=="Between",Dice,NA)
  MDiceWithin <- mean(Dice.Within[lower.tri(Dice.Within, diag = FALSE)],na.rm=TRUE)
  MDiceBetween <- mean(Dice.Between[lower.tri(Dice.Between, diag = FALSE)],na.rm=TRUE)
  dat.Perm[k] <- MDiceWithin - MDiceBetween 
  flush.console()
  if(k %% 10 == 0) { print(paste0("Finished ", k, " out of ", nPerm, " simulations")) } 
}

hist(dat.Perm)
abline(v=Emp.MDiceWithin - Emp.MDiceBetween, col="red") 
pSRG <- (sum(abs(dat.Perm) >= abs(Emp.MDiceWithin - Emp.MDiceBetween)) + 1) / (nPerm + 1)
# P-value = 0.000999


# Significance thresholds are P>=0.975 and P<=0.025. 
# This is because we are looking at the deviation from 0, which can either be negative or positive.
# Because the distribution of differences is not necessarily symmetric around zero, we cannot calculate P-values using absolute values.


# Repertoire similarity within group A (A-A), within group B (B-B), within group C2, within group C4, within group L, within group K

# step 2b: define within and between setting diads --- for A, B, C2, C4, L and K groups 
## create data frame with Subject as rows and columns
dat.wth.btw2 <- data.frame(matrix(rep(NA, nrow(IndRepSize) * nrow(IndRepSize)), ncol=nrow(IndRepSize)))
colnames(dat.wth.btw2) <- IndRepSize$Subject
rownames(dat.wth.btw2) <- IndRepSize$Subject
## loop over all pairs of individuals
for(i in 1:(nrow(IndRepSize))) {
  for(j in 1:(nrow(IndRepSize))) {
    wth.btw2 <- NA
    Ind1 <- subset(Data, Subject==IndRepSize$Subject[i])
    Ind2 <- subset(Data, Subject==IndRepSize$Subject[j])
    if(length(table(Ind1$Group))>1) { print("Ohoh, better check my data again!") }
    if(length(table(Ind2$Group))>1) { print("Ohoh, better check my data again!") }
    if(Ind1$Group[1]=="A" & Ind2$Group[1]=="A") { wth.btw2 <- "A" }
    if(Ind1$Group[1]=="B" & Ind2$Group[1]=="B") { wth.btw2 <- "B" }
    if(Ind1$Group[1]=="C2" & Ind2$Group[1]=="C2") { wth.btw2 <- "C2" }
    if(Ind1$Group[1]=="C4" & Ind2$Group[1]=="C4") { wth.btw2 <- "C4" }
    if(Ind1$Group[1]=="L" & Ind2$Group[1]=="L") { wth.btw2 <- "L" }
    if(Ind1$Group[1]=="K" & Ind2$Group[1]=="K") { wth.btw2 <- "K" }
    dat.wth.btw2[which(rownames(dat.wth.btw2)==IndRepSize$Subject[i]),which(colnames(dat.wth.btw2)==IndRepSize$Subject[j])] <- wth.btw2
  }
}
## set diagonal to NA
wth.btw2 <- as.matrix(dat.wth.btw2)
diag(wth.btw2) <- NA
#wth.btw2 <- as.data.frame(wth.btw2)

# step 3b: select within and between group diads from the Dice data frame ---
Dice.A <- ifelse(wth.btw2=="A",Dice,NA)
Dice.B <- ifelse(wth.btw2=="B",Dice,NA)
Dice.C2 <- ifelse(wth.btw2=="C2",Dice,NA)
Dice.C4 <- ifelse(wth.btw2=="C4",Dice,NA)
Dice.L <- ifelse(wth.btw2=="L",Dice,NA)
Dice.K <- ifelse(wth.btw2=="K",Dice,NA)

(Emp.MDiceA <- mean(Dice.A[lower.tri(Dice.A, diag = FALSE)],na.rm=TRUE))
# Mean Dice within group A = 0.5968785 --> big individual variation (0 =  no overlap ; 1 = 100% overlap)
(Emp.MDiceB <- mean(Dice.B[lower.tri(Dice.B, diag = FALSE)],na.rm=TRUE))
# Mean Dice within group B = 0.4133377
(Emp.MDiceC2 <- mean(Dice.C2[lower.tri(Dice.C2, diag = FALSE)],na.rm=TRUE))
# Mean Dice within group C2 = 0.5885302
(Emp.MDiceC4 <- mean(Dice.C4[lower.tri(Dice.C4, diag = FALSE)],na.rm=TRUE))
# Mean Dice within group C4 = 0.5648174
(Emp.MDiceL <- mean(Dice.L[lower.tri(Dice.L, diag = FALSE)],na.rm=TRUE))
# Mean Dice within group L = 0.4897272
(Emp.MDiceK <- mean(Dice.K[lower.tri(Dice.K, diag = FALSE)],na.rm=TRUE))
# Mean Dice within group K = 0.561168

## step 4b: Combiner toutes les infos dans un tableau (1 ligne = 1 diade --> create data frame dd with columns Ind1, Ind2, Dice, A/B)
## Tranform A matrix in table
yA <- expand.grid(rownames(Dice.A), colnames(Dice.A)) 
labsA <- yA[as.vector(upper.tri(Dice.A, diag = FALSE)), ] 
yA <- cbind(labsA, Dice.A[upper.tri(Dice.A,diag=FALSE)]) 
colnames(yA) <- c("Ind1","Ind2","Dice")
yA <- yA[!is.na(yA$Dice), ] 
yA$Group <- rep("A",nrow(yA)) 

## Tranform B matrix in table
yB <- expand.grid(rownames(Dice.B), colnames(Dice.B))
labsB <- yB[as.vector(upper.tri(Dice.B, diag = FALSE)), ]
yB <- cbind(labsB, Dice.B[upper.tri(Dice.B,diag=FALSE)])
colnames(yB) <- c("Ind1","Ind2","Dice")
yB <- yB[!is.na(yB$Dice), ]
yB$Group <- rep("B",nrow(yB))

## Tranform C2 matrix in table
yC2 <- expand.grid(rownames(Dice.C2), colnames(Dice.C2))
labsC2 <- yC2[as.vector(upper.tri(Dice.C2, diag = FALSE)), ] 
yC2 <- cbind(labsC2, Dice.C2[upper.tri(Dice.C2,diag=FALSE)]) 
colnames(yC2) <- c("Ind1","Ind2","Dice")
yC2 <- yC2[!is.na(yC2$Dice), ] 
yC2$Group <- rep("C2",nrow(yC2)) 

## Tranform C4 matrix in table
yC4 <- expand.grid(rownames(Dice.C4), colnames(Dice.C4)) 
labsC4 <- yC4[as.vector(upper.tri(Dice.C4, diag = FALSE)), ] 
yC4 <- cbind(labsC4, Dice.C4[upper.tri(Dice.C4,diag=FALSE)]) 
colnames(yC4) <- c("Ind1","Ind2","Dice")
yC4 <- yC4[!is.na(yC4$Dice), ] 
yC4$Group <- rep("C4",nrow(yC4)) 

## Tranform L matrix in table
yL <- expand.grid(rownames(Dice.L), colnames(Dice.L)) 
labsL <- yL[as.vector(upper.tri(Dice.L, diag = FALSE)), ] 
yL <- cbind(labsL, Dice.L[upper.tri(Dice.L,diag=FALSE)]) 
colnames(yL) <- c("Ind1","Ind2","Dice")
yL <- yL[!is.na(yL$Dice), ] 
yL$Group <- rep("L",nrow(yL)) 

## Tranform K matrix in table
yK <- expand.grid(rownames(Dice.K), colnames(Dice.K)) 
labsK <- yK[as.vector(upper.tri(Dice.K, diag = FALSE)), ] 
yK <- cbind(labsK, Dice.K[upper.tri(Dice.K,diag=FALSE)]) 
colnames(yK) <- c("Ind1","Ind2","Dice")
yK <- yK[!is.na(yK$Dice), ] 
yK$Group <- rep("K",nrow(yK)) 

## Combine A/B/C2/C4/L/K tables
dd <- rbind(yA,yB,yC2,yC4,yL,yK)

Data.sub <- Data[,c("Subject","Group")] 
Data.sub <- Data.sub[!duplicated(Data.sub), ] 
dd <- merge(dd, Data.sub, by.x="Ind1", by.y="Subject", sort=FALSE, all.y=FALSE) 
dd <- merge(dd, Data.sub, by.x="Ind2", by.y="Subject", sort=FALSE) 
colnames(dd)[which(colnames(dd)=="Group.x")] <- "Group.ind1"
colnames(dd)[which(colnames(dd)=="Group.y")] <- "Group.ind2"
# Order columns
dd<-dd[,c("Ind1","Ind2","Dice","Group","Group.ind1","Group.ind2")]

# step 5b: matrix permutation test ---
nPerm <- 1000
dat.Perm2 <- c()
for(k in 1:nPerm) {
  ## permute group in Data
  ind.group.perm <- Data.sub
  ind.group.perm$Group <- sample(ind.group.perm$Group, nrow(ind.group.perm), replace=FALSE)
  Data.inc.perm <- Data[,-which(colnames(Data)=="Group")] 
  Data.inc.perm <- merge(Data.inc.perm, ind.group.perm, by.x="Subject", by.y="Subject", sort=FALSE)
  ## repeat step 2 & 3
  ## create data frame with Subject as rows and columns (step2)
  dat.wth.btw2 <- data.frame(matrix(rep(NA, nrow(IndRepSize) * nrow(IndRepSize)), ncol=nrow(IndRepSize)))
  colnames(dat.wth.btw2) <- IndRepSize$Subject
  rownames(dat.wth.btw2) <- IndRepSize$Subject
  ## loop over all pairs of individuals (step2)
  for(i in 1:(nrow(IndRepSize))) {
    for(j in 1:(nrow(IndRepSize))) {
      Ind1 <- subset(Data.inc.perm, Subject==IndRepSize$Subject[i])
      Ind2 <- subset(Data.inc.perm, Subject==IndRepSize$Subject[j])
      if(length(table(Ind1$Group))>1) { print("Ohoh, better check my data again!") }
      if(length(table(Ind2$Group))>1) { print("Ohoh, better check my data again!") }
      if(Ind1$Group[1]=="A" & Ind2$Group[1]=="A") { wth.btw2 <- "A" }
      if(Ind1$Group[1]=="B" & Ind2$Group[1]=="B") { wth.btw2 <- "B" }
      if(Ind1$Group[1]=="C2" & Ind2$Group[1]=="C2") { wth.btw2 <- "C2" }
      if(Ind1$Group[1]=="C4" & Ind2$Group[1]=="C4") { wth.btw2 <- "C4" }
      if(Ind1$Group[1]=="L" & Ind2$Group[1]=="L") { wth.btw2 <- "L" }
      if(Ind1$Group[1]=="K" & Ind2$Group[1]=="K") { wth.btw2 <- "K" }
      dat.wth.btw2[which(rownames(dat.wth.btw2)==IndRepSize$Subject[i]),which(colnames(dat.wth.btw2)==IndRepSize$Subject[j])] <- wth.btw2
    }
  }
  ## set diagonal to NA (step6)
  wth.btw2 <- as.matrix(dat.wth.btw2)
  diag(wth.btw2) <- NA
  Dice.A <- ifelse(wth.btw2=="A",Dice,NA)
  Dice.B <- ifelse(wth.btw2=="B",Dice,NA)
  Dice.C2 <- ifelse(wth.btw2=="C2",Dice,NA)
  Dice.C4 <- ifelse(wth.btw2=="C4",Dice,NA)
  Dice.L <- ifelse(wth.btw2=="L",Dice,NA)
  Dice.K <- ifelse(wth.btw2=="K",Dice,NA)
  MDiceA <- mean(Dice.A[lower.tri(Dice.A, diag = FALSE)],na.rm=TRUE)
  MDiceB <- mean(Dice.B[lower.tri(Dice.B, diag = FALSE)],na.rm=TRUE)
  MDiceC2 <- mean(Dice.C2[lower.tri(Dice.C2, diag = FALSE)],na.rm=TRUE)
  MDiceC4 <- mean(Dice.C4[lower.tri(Dice.C4, diag = FALSE)],na.rm=TRUE)
  MDiceL <- mean(Dice.L[lower.tri(Dice.L, diag = FALSE)],na.rm=TRUE)
  MDiceK <- mean(Dice.K[lower.tri(Dice.K, diag = FALSE)],na.rm=TRUE)
  # omnibus statistic
  mu <- c(MDiceA, MDiceB, MDiceC2, MDiceC4, MDiceL, MDiceK)
  dat.Perm2[k] <- sum((mu - mean(mu))^2)
  # pairwise differences
  perm.pairwise[k, ] <- c(
    MDiceA - MDiceB,
    MDiceA - MDiceC2,
    MDiceA - MDiceC4,
    MDiceA - MDiceL,
    MDiceA - MDiceK,
    MDiceB - MDiceC2,
    MDiceB - MDiceC4,
    MDiceB - MDiceL,
    MDiceB - MDiceK,
    MDiceC2 - MDiceC4,
    MDiceC2 - MDiceL,
    MDiceC2 - MDiceK,
    MDiceC4 - MDiceL,
    MDiceC4 - MDiceK,
    MDiceL - MDiceK
  )
  flush.console()
  if(k %% 10 == 0) { print(paste0("Finished ", k, " out of ", nPerm, " simulations")) } 
}
hist(dat.Perm2)
abline(v=Emp.MDiceA - Emp.MDiceB - Emp.MDiceC2 - Emp.MDiceC4 - Emp.MDiceL - Emp.MDiceK, col="red") 

mu.emp <- c(Emp.MDiceA, Emp.MDiceB, Emp.MDiceC2, Emp.MDiceC4, Emp.MDiceL, Emp.MDiceK)
mean.emp <- mean(mu.emp)
SS.emp <- sum((mu.emp - mean.emp)^2) 

p.valueSRG <- (sum(dat.Perm2 >= SS.emp) + 1) / (length(dat.Perm2) + 1) #P-value = 0.000999

#only if the variance differs:
pairwise.emp <- c(
  Emp.MDiceA - Emp.MDiceB,
  Emp.MDiceA - Emp.MDiceC2,
  Emp.MDiceA - Emp.MDiceC4,
  Emp.MDiceA - Emp.MDiceL,
  Emp.MDiceA - Emp.MDiceK,
  Emp.MDiceB - Emp.MDiceC2,
  Emp.MDiceB - Emp.MDiceC4,
  Emp.MDiceB - Emp.MDiceL,
  Emp.MDiceB - Emp.MDiceK,
  Emp.MDiceC2 - Emp.MDiceC4,
  Emp.MDiceC2 - Emp.MDiceL,
  Emp.MDiceC2 - Emp.MDiceK,
  Emp.MDiceC4 - Emp.MDiceL,
  Emp.MDiceC4 - Emp.MDiceK,
  Emp.MDiceL - Emp.MDiceK
)

pvals <- sapply(1:15, function(i) {
  valid <- !is.na(perm.pairwise[, i])
  (sum(abs(perm.pairwise[valid, i]) >= abs(pairwise.emp[i])) + 1) /
    (sum(valid) + 1)
})

pvals.holm <- p.adjust(pvals, method = "holm")

contrast.names <- c(
  "A-B","A-C2","A-C4","A-L","A-K",
  "B-C2","B-C4","B-L","B-K",
  "C2-C4","C2-L","C2-K",
  "C4-L","C4-K",
  "L-K"
)

results <- data.frame(
  Contrast = contrast.names,
  p_raw = pvals,
  p_holm = pvals.holm
)

results
#   Contrast       p_raw     p_holm
#1       A-B 0.000999001 0.01498501 *
#2      A-C2 0.770229770 1.00000000
#3      A-C4 0.362637363 1.00000000
#4       A-L 0.000999001 0.01498501 *
#5       A-K 0.187812188 1.00000000
#6      B-C2 0.000999001 0.01498501 *
#7      B-C4 0.001998002 0.01798202 *
#8       B-L 0.050949051 0.35664336
#9       B-K 0.000999001 0.01498501 *
#10    C2-C4 0.465534466 1.00000000
#11     C2-L 0.000999001 0.01498501 *
#12     C2-K 0.177822178 1.00000000
#13     C4-L 0.011988012 0.09590410
#14     C4-K 0.895104895 1.00000000
#15      L-K 0.000999001 0.01498501 *

# Boxplot Dice within/between (steps a)

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

levels(WBdd$Group.ind1) <- c("A", "B", "L", "C2", "C4", "K") # dice = dd ?

levels(dd$Group.ind1) <- c("A", "B", "L", "C2", "C4", "K") # dice = dd ?

WBddWithin <- subset(WBdd, WBdd$WithinBetween=="Within")


F1 <- ggplot() + 
  geom_boxplot(WBddWithin, mapping=aes(x = Group.ind1, y = Dice), width = 0.9, fill = c("#007ABB", "#00AFBB","#259C39", "#E7B800","#E79A00", "#ba0f09")) +
  geom_boxplot(WBdd, mapping=aes(x = WithinBetween, y = Dice), width = 0.9, fill = c("grey90","grey90")) +
  geom_point(WBddWithin, mapping=aes(x = Group.ind1, y = Dice), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_point(WBdd, mapping=aes(x = WithinBetween, y = Dice), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "grey40") + 
  theme_angele_ss +
  ggtitle("(a)") +
  scale_y_continuous("Repertoire similarity among individuals", breaks = seq(0, 1, by = 0.2)) +
  scale_x_discrete(" ",
                   limits = c("A", "B", "L", "C2", "C4", "K", "Between", "Within"),
                   labels = c("Within \ngroup A", "Within \ngroup B", "Within \ngroup L", "Within \ngroup C2", "Within \ngroup C4", "Within \ngroup K", "Between \ngroups", "Within \ngroups"))+
  #facet_wrap(~Dice, scales='free_x')+
  stat_summary(WBddWithin, mapping=aes(x = Group.ind1, y = Dice), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3) +
  stat_summary(WBdd, mapping=aes(x = WithinBetween, y = Dice), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3)

# + theme(axis.text.x = element_blank()) 

print(F1)
