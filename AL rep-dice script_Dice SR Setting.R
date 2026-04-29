5. SETTING SIGNAL REPERTOIRE SIMILARITY 

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

# Omit all levels of Subject that contributed fewer than 30 cases --- 118 --> 89 individuals
nb.signals.per.ind <- data.frame(Data %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(nb.signals.per.ind) <- c("Subject","NsignalsTot") 
nb.signals.per.ind$Subject <- as.character(nb.signals.per.ind$Subject)
nb.signals.per.ind <- subset(nb.signals.per.ind, nb.signals.per.ind$NsignalsTot>=30) # Keeps only the values in "column" that appear more than 30 times (N >= 30)
Data <- subset(Data, Subject %in% nb.signals.per.ind$Subject)

Data <- Data[!is.na(Data$Group), ]

IndRepSize <- aggregate(data=Data, Behavior ~ Subject, function(x) length(unique(x))) ; colnames(IndRepSize) <- c("Subject","RepertoireSize") 
dat.ind.bhv2 <- data.frame(Data %>% group_by(Subject,Behavior) %>% dplyr::summarize(count=n())) ; colnames(dat.ind.bhv2) <- c("Subject","Behavior","NsignalsPbehavior") 

________________________________________________________________________________________________________________________________________________________________________________________________


# Individual repertoire similarity (Captive-Captive/Sanctuary-Sanctuary, Captive-Sanctuary): 

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

# step 2: define within- and between-setting dyads 
## create data frame with Subject as rows and columns
SETdat.wth.btw <- data.frame(matrix(rep(NA, nrow(IndRepSize) * nrow(IndRepSize)), ncol=nrow(IndRepSize)))
colnames(SETdat.wth.btw) <- IndRepSize$Subject
rownames(SETdat.wth.btw) <- IndRepSize$Subject
## loop over all pairs of individuals
for(i in 1:(nrow(IndRepSize))) {
  for(j in 1:(nrow(IndRepSize))) {
    Ind1 <- subset(Data, Subject==IndRepSize$Subject[i])
    Ind2 <- subset(Data, Subject==IndRepSize$Subject[j])
    if(length(table(Ind1$Setting))>1) { print("Ohoh, better check my data again!") }
    if(length(table(Ind2$Setting))>1) { print("Ohoh, better check my data again!") }
    if(Ind1$Setting[1]==Ind2$Setting[1]) { SETwth.btw <- "Within" } else { SETwth.btw <- "Between" }
    SETdat.wth.btw[which(rownames(SETdat.wth.btw)==IndRepSize$Subject[i]),which(colnames(SETdat.wth.btw)==IndRepSize$Subject[j])] <- SETwth.btw
  }
}
## set diagonal to NA
SETwth.btw <- as.matrix(SETdat.wth.btw)
diag(SETwth.btw) <- NA
#SETwth.btw <- as.data.frame(SETwth.btw)

# step 3: select within and between setting dyads from the Dice data frame (N=19)
Dice.WithinSET <- ifelse(SETwth.btw=="Within",Dice,NA)
Dice.BetweenSET <- ifelse(SETwth.btw=="Between",Dice,NA)

(Emp.MDiceWithinSET <- mean(Dice.WithinSET[lower.tri(Dice.WithinSET, diag = FALSE)],na.rm=TRUE)) 
# Mean Dice within settings (settings C, S et W confondus) = 0.5232495
(Emp.MDiceBetweenSET <- mean(Dice.BetweenSET[lower.tri(Dice.BetweenSET, diag = FALSE)],na.rm=TRUE)) 
# Mean Dice between (C-S) settings = 0.4665083

# Repertoire similarity between groups Within et Between

# step 4
## Tranform Within matrix in table
yWithin <- expand.grid(rownames(Dice.WithinSET), colnames(Dice.WithinSET)) 
labsWithin <- yWithin[as.vector(upper.tri(Dice.WithinSET, diag = FALSE)), ] 
yWithin <- cbind(labsWithin, Dice.WithinSET[upper.tri(Dice.WithinSET,diag=FALSE)]) 
colnames(yWithin) <- c("Ind1","Ind2","Dice")
yWithin <- yWithin[!is.na(yWithin$Dice), ] 
yWithin$WithinBetween <- rep("WithinS",nrow(yWithin)) 
mean(yWithin$Dice)

## Tranform Between matrix in table
yBetween <- expand.grid(rownames(Dice.BetweenSET), colnames(Dice.BetweenSET))
labsBetween <- yBetween[as.vector(upper.tri(Dice.BetweenSET, diag = FALSE)), ]
yBetween <- cbind(labsBetween, Dice.BetweenSET[upper.tri(Dice.BetweenSET,diag=FALSE)])
colnames(yBetween) <- c("Ind1","Ind2","Dice")
yBetween <- yBetween[!is.na(yBetween$Dice), ]
yBetween$WithinBetween <- rep("BetweenS",nrow(yBetween))
mean(yBetween$Dice)

##Combine Within and Between tables
WBddSET <- rbind(yWithin,yBetween)

Data.sub <- Data[,c("Subject","Setting")] 
Data.sub <- Data.sub[!duplicated(Data.sub), ] 
WBddSET <- merge(WBddSET, Data.sub, by.x="Ind1", by.y="Subject", sort=FALSE, all.y=FALSE) 
WBddSET <- merge(WBddSET, Data.sub, by.x="Ind2", by.y="Subject", sort=FALSE) 
colnames(WBddSET)[which(colnames(WBddSET)=="Setting.x")] <- "Setting.ind1"
colnames(WBddSET)[which(colnames(WBddSET)=="Setting.y")] <- "Setting.ind2"
# V?rifier l'ordre des colonnes
WBddSET<-WBddSET[,c("Ind1","Ind2","Dice","WithinBetween","Setting.ind1","Setting.ind2")]

# step 5: matrix permutation test ---
nPerm <- 1000
dat.Perm <- c()
for(k in 1:nPerm) {
  ## permute Setting in Data
  ind.Setting.perm <- Data.sub
  ind.Setting.perm$Setting <- sample(ind.Setting.perm$Setting, nrow(ind.Setting.perm), replace=FALSE)
  Data.inc.perm <- Data[,-which(colnames(Data)=="Setting")] 
  Data.inc.perm <- merge(Data.inc.perm, ind.Setting.perm, by.x="Subject", by.y="Subject", sort=FALSE)
  ## repeat step 2 & 3
  ## create data frame with animal_id as rows and columns (step2)
  SETdat.wth.btw <- data.frame(matrix(rep(NA, nrow(IndRepSize) * nrow(IndRepSize)), ncol=nrow(IndRepSize)))
  colnames(SETdat.wth.btw) <- IndRepSize$Subject
  rownames(SETdat.wth.btw) <- IndRepSize$Subject
  ## loop over all pairs of individuals (step2)
  for(i in 1:(nrow(IndRepSize))) {
    for(j in 1:(nrow(IndRepSize))) {
      Ind1 <- subset(Data.inc.perm, Subject==IndRepSize$Subject[i])
      Ind2 <- subset(Data.inc.perm, Subject==IndRepSize$Subject[j])
      if(length(table(Ind1$Setting))>1) { print("Ohoh, better check my data again!") }
      if(length(table(Ind2$Setting))>1) { print("Ohoh, better check my data again!") }
      if(Ind1$Setting[1]==Ind2$Setting[1]) { SETwth.btw <- "Within" } else { SETwth.btw <- "Between" }
      SETdat.wth.btw[which(rownames(SETdat.wth.btw)==IndRepSize$Subject[i]),which(colnames(SETdat.wth.btw)==IndRepSize$Subject[j])] <- SETwth.btw
    }
  }
  ## set diagonal to NA (step3)
  SETwth.btw <- as.matrix(SETdat.wth.btw)
  diag(SETwth.btw) <- NA
  Dice.WithinSET <- ifelse(SETwth.btw=="Within",Dice,NA)
  Dice.BetweenSET <- ifelse(SETwth.btw=="Between",Dice,NA)
  MDiceWithinSET <- mean(Dice.WithinSET[lower.tri(Dice.WithinSET, diag = FALSE)],na.rm=TRUE)
  MDiceBetweenSET <- mean(Dice.BetweenSET[lower.tri(Dice.BetweenSET, diag = FALSE)],na.rm=TRUE)
  dat.Perm[k] <- MDiceWithinSET - MDiceBetweenSET 
  flush.console()
  if(k %% 10 == 0) { print(paste0("Finished ", k, " out of ", nPerm, " simulations")) } 
}
hist(dat.Perm)
abline(v=Emp.MDiceWithinSET - Emp.MDiceBetweenSET, col="red") 
pSRS <- (sum(abs(dat.Perm) >= abs(Emp.MDiceWithinSET - Emp.MDiceBetweenSET)) + 1) / (nPerm + 1)
# P-value = 0.000999


# Significance thresholds are P>=0.975 and P<=0.025. 
# This is because we are looking at the deviation from 0, which can either be negative or positive.
# Because the distribution of differences is not necessarily symmetric around zero, we cannot calculate P-values using absolute values.


# Repertoire similarity within setting captive (C-C), within setting sanctuary (S-S) & within setting wild (W-W)

# step 2b: define within and between setting diads --- for C and S settings 
## create data frame with animal_id as rows and columns
SETdat.wth.btw2 <- data.frame(matrix(rep(NA, nrow(IndRepSize) * nrow(IndRepSize)), ncol=nrow(IndRepSize)))
colnames(SETdat.wth.btw2) <- IndRepSize$Subject
rownames(SETdat.wth.btw2) <- IndRepSize$Subject
## loop over all pairs of individuals
for(i in 1:(nrow(IndRepSize))) {
  for(j in 1:(nrow(IndRepSize))) {
    SETwth.btw2 <- NA
    Ind1 <- subset(Data, Subject==IndRepSize$Subject[i])
    Ind2 <- subset(Data, Subject==IndRepSize$Subject[j])
    if(length(table(Ind1$Setting))>1) { print("Ohoh, better check my data again!") }
    if(length(table(Ind2$Setting))>1) { print("Ohoh, better check my data again!") }
    if(Ind1$Setting[1]=="Captive" & Ind2$Setting[1]=="Captive") { SETwth.btw2 <- "Captive" }
    if(Ind1$Setting[1]=="Sanctuary" & Ind2$Setting[1]=="Sanctuary") { SETwth.btw2 <- "Sanctuary" }
    if(Ind1$Setting[1]=="Wild" & Ind2$Setting[1]=="Wild") { SETwth.btw2 <- "Wild" }
    SETdat.wth.btw2[which(rownames(SETdat.wth.btw2)==IndRepSize$Subject[i]),which(colnames(SETdat.wth.btw2)==IndRepSize$Subject[j])] <- SETwth.btw2
  }
}
## set diagonal to NA
SETwth.btw2 <- as.matrix(SETdat.wth.btw2)
diag(SETwth.btw2) <- NA
#SETwth.btw2 <- as.data.frame(SETwth.btw2)

# step 3b: select within and between group diads from the Dice data frame ---
Dice.C <- ifelse(SETwth.btw2=="Captive",Dice,NA)
Dice.S <- ifelse(SETwth.btw2=="Sanctuary",Dice,NA)
Dice.W <- ifelse(SETwth.btw2=="Wild",Dice,NA)

(Emp.MDiceC <- mean(Dice.C[lower.tri(Dice.C, diag = FALSE)],na.rm=TRUE))
# Mean Dice within Captive = 0.4922726 --> big individual variation (0 =  no overlap ; 1 = 100% overlap)
(Emp.MDiceS <- mean(Dice.S[lower.tri(Dice.S, diag = FALSE)],na.rm=TRUE))
# Mean Dice within Sanctuary = 0.5513234
(Emp.MDiceW <- mean(Dice.W[lower.tri(Dice.W, diag = FALSE)],na.rm=TRUE))
# Mean Dice within Wild = 0.561168

## step 4b: 
## Tranform C matrix in table
yC <- expand.grid(rownames(Dice.C), colnames(Dice.C)) 
labsC <- yC[as.vector(upper.tri(Dice.C, diag = FALSE)), ] 
yC <- cbind(labsC, Dice.C[upper.tri(Dice.C,diag=FALSE)]) 
colnames(yC) <- c("Ind1","Ind2","Dice")
yC <- yC[!is.na(yC$Dice), ] 
yC$SET <- rep("Captive",nrow(yC)) 
mean(yC$Dice)

## Tranform S matrix in table
yS <- expand.grid(rownames(Dice.S), colnames(Dice.S))
labsS <- yS[as.vector(upper.tri(Dice.S, diag = FALSE)), ]
yS <- cbind(labsS, Dice.S[upper.tri(Dice.S,diag=FALSE)])
colnames(yS) <- c("Ind1","Ind2","Dice")
yS <- yS[!is.na(yS$Dice), ]
yS$SET <- rep("Sanctuary",nrow(yS))
mean(yS$Dice)

## Tranform W matrix in table
yW <- expand.grid(rownames(Dice.W), colnames(Dice.W))
labsW <- yW[as.vector(upper.tri(Dice.W, diag = FALSE)), ]
yW <- cbind(labsW, Dice.W[upper.tri(Dice.W,diag=FALSE)])
colnames(yW) <- c("Ind1","Ind2","Dice")
yW <- yW[!is.na(yW$Dice), ]
yW$SET <- rep("Wild",nrow(yW))
mean(yW$Dice)

## Combine C, S and W tables
SETdd <- rbind(yC,yS,yW)

Data.sub <- Data[,c("Subject","Setting")] 
Data.sub <- Data.sub[!duplicated(Data.sub), ] 
SETdd <- merge(SETdd, Data.sub, by.x="Ind1", by.y="Subject", sort=FALSE, all.y=FALSE) 
SETdd <- merge(SETdd, Data.sub, by.x="Ind2", by.y="Subject", sort=FALSE) 
colnames(SETdd)[which(colnames(SETdd)=="Setting.x")] <- "Setting.ind1"
colnames(SETdd)[which(colnames(SETdd)=="Setting.y")] <- "Setting.ind2"
# V?rifier l'ordre des colonnes
SETdd<-SETdd[,c("Ind1","Ind2","Dice","SET","Setting.ind1","Setting.ind2")]

mu.empSRS <- c(Emp.MDiceC, Emp.MDiceS, Emp.MDiceW)
names(mu.empSRS) <- c("Captive","Sanctuary","Wild")
# Get the number of individuals per group
group_sizesSET <- table(Data.sub$Setting)[c("Captive","Sanctuary","Wild")]
wSET <- group_sizesSET * (group_sizesSET - 1) / 2
# Compute the weighted sums of squares
mean.w.empSRS <- weighted.mean(mu.empSRS, wSET=wSET, na.rm=TRUE)
SS.empSRS <- sum(wSET * (mu.empSRS - mean.w.empSRS)^2, na.rm=TRUE) 


# step 5b: matrix permutation test ---
nPerm <- 1000
dat.Perm3 <- c()
perm.pairwiseSRS <- matrix(NA_real_, nrow = nPerm, ncol = 3)
for(k in 1:nPerm) {
  ## permute group in Data
  ind.Setting.perm <- Data.sub
  ind.Setting.perm$Setting <- sample(ind.Setting.perm$Setting, nrow(ind.Setting.perm), replace=FALSE)
  Data.inc.perm <- Data[,-which(colnames(Data)=="Setting")] 
  Data.inc.perm <- merge(Data.inc.perm, ind.Setting.perm, by.x="Subject", by.y="Subject", sort=FALSE)
  ## repeat step 2 & 3
  ## create data frame with animal_id as rows and columns (step2)
  SETdat.wth.btw2 <- data.frame(matrix(rep(NA, nrow(IndRepSize) * nrow(IndRepSize)), ncol=nrow(IndRepSize)))
  colnames(SETdat.wth.btw2) <- IndRepSize$Subject
  rownames(SETdat.wth.btw2) <- IndRepSize$Subject
  ## loop over all pairs of individuals (step2)
  for(i in 1:(nrow(IndRepSize))) {
    for(j in 1:(nrow(IndRepSize))) {
      Ind1 <- subset(Data.inc.perm, Subject==IndRepSize$Subject[i])
      Ind2 <- subset(Data.inc.perm, Subject==IndRepSize$Subject[j])
      if(length(table(Ind1$Setting))>1) { print("Ohoh, better check my data again!") }
      if(length(table(Ind2$Setting))>1) { print("Ohoh, better check my data again!") }
      if(Ind1$Setting[1]=="Captive" & Ind2$Setting[1]=="Captive") { SETwth.btw2 <- "Captive" }
      if(Ind1$Setting[1]=="Sanctuary" & Ind2$Setting[1]=="Sanctuary") { SETwth.btw2 <- "Sanctuary" }
      if(Ind1$Setting[1]=="Wild" & Ind2$Setting[1]=="Wild") { SETwth.btw2 <- "Wild" }
      SETdat.wth.btw2[which(rownames(SETdat.wth.btw2)==IndRepSize$Subject[i]),which(colnames(SETdat.wth.btw2)==IndRepSize$Subject[j])] <- SETwth.btw2
    }
  }
  ## set diagonal to NA (step6)
  SETwth.btw2 <- as.matrix(SETdat.wth.btw2)
  diag(SETwth.btw2) <- NA
  Dice.C <- ifelse(SETwth.btw2=="Captive",Dice,NA)
  Dice.S <- ifelse(SETwth.btw2=="Sanctuary",Dice,NA)
  Dice.W <- ifelse(SETwth.btw2=="Wild",Dice,NA)
  MDiceC <- mean(Dice.C[lower.tri(Dice.C, diag = FALSE)],na.rm=TRUE)
  MDiceS <- mean(Dice.S[lower.tri(Dice.S, diag = FALSE)],na.rm=TRUE)
  MDiceW <- mean(Dice.W[lower.tri(Dice.W, diag = FALSE)],na.rm=TRUE)
  # omnibus statistic
  muSRG <- c(MDiceC, MDiceS, MDiceW)
  mean.wSRG <- weighted.mean(muSRG, wSET=wSET, na.rm=TRUE)
  dat.Perm3[k] <- sum(wSET * (muSRG - mean.wSRG)^2, na.rm=TRUE)
  # pairwise differences
  perm.pairwiseSRS[k, ] <- c(
    MDiceC - MDiceS,
    MDiceC - MDiceW,
    MDiceS - MDiceW
  )
  flush.console()
  if(k %% 10 == 0) { print(paste0("Finished ", k, " out of ", nPerm, " simulations")) } 
}
hist(dat.Perm3)
abline(v=SS.empSRS, col="red") 

p.valueSRS <- (sum(dat.Perm3 >= SS.empSRS) + 1) / (length(dat.Perm3) + 1) #P-value = 0.000999

#only if the variance differs:
pairwiseSRS.emp <- c(
  Emp.MDiceC - Emp.MDiceS,
  Emp.MDiceC - Emp.MDiceW,
  Emp.MDiceS - Emp.MDiceW
)

pvalsSRS <- sapply(1:3, function(i) {
  valid <- !is.na(perm.pairwiseSRS[, i])
  (sum(abs(perm.pairwiseSRS[valid, i]) >= abs(pairwiseSRS.emp[i])) + 1) /
    (sum(valid) + 1)
})

pvalsSRS.holm <- p.adjust(pvalsSRS, method = "holm")

contrast.namesSET <- c(
  "C-S","C-W","S-W")

resultsSRS <- data.frame(
  Contrast = contrast.namesSET,
  p_raw = pvalsSRS,
  p_holm = pvalsSRS.holm
)

resultsSRS
#  Contrast       p_raw      p_holm
#1      C-S 0.000999001 0.002997003 *
#2      C-W 0.001998002 0.003996004 *
#3      S-W 0.625374625 0.625374625 


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

levels(WBddSET$Setting.ind1) <- c("Captive", "Sanctuary", "Wild") # Dice = dd ?

WBddSETWithin <- subset(WBddSET, WBddSET$WithinBetween=="WithinS")

F2 <- ggplot() + 
  geom_boxplot(WBddSETWithin, mapping=aes(x = Setting.ind1, y = Dice), width = 0.9, fill = c("#068591", "#ba8211", "#ba0f09")) +
  geom_boxplot(WBddSET, mapping=aes(x = WithinBetween, y = Dice), width = 0.9, fill = c("grey90","grey90")) +
  geom_point(WBddSETWithin, mapping=aes(x = Setting.ind1, y = Dice), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_point(WBddSET, mapping=aes(x = WithinBetween, y = Dice), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey40") + 
  theme_angele_ss +
  ggtitle("(b)") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_x_discrete(" ",
                   limits = c("Captive", "Sanctuary", "Wild", "BetweenS", "WithinS"),
                   labels = c("Within \nZoo", "Within \nSanctuary", "Within \nWild", "Between \nsettings", "Within \nsettings"))+
  #facet_wrap(~Dice, scales='free_x')+
  stat_summary(WBddSETWithin, mapping=aes(x = Setting.ind1, y = Dice), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3) +
  stat_summary(WBddSET, mapping=aes(x = WithinBetween, y = Dice), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

print(F2)
ggarrange(F1, F2, widths = c(1.7,1))

