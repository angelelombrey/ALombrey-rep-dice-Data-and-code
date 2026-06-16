## 4. GROUP MULTICOMPONENT REPERTOIRE SIMILARITY ##

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
path <- "D:\\Documents\\CoAuthoredManuscripts\\2026__Lombrey__Chimpanzee_Plasticity\\DataAnalysis\\"
library(openxlsx)
Data <- read.xlsx(xlsxFile=paste0(path,"Data.xlsx"), sheet=1, colNames=TRUE)

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
SamplingEffort <- data.frame(datMultiRep %>% group_by(Subject) %>% dplyr::summarize(count=n(), .groups = "drop")) ; colnames(SamplingEffort) <- c("Subject","TotComActs") 
datMultiRep <- merge(SamplingEffort, datMultiRep, by=c("Subject"))

datMultiRep <- datMultiRep[!is.na(datMultiRep$Group), ]

datMultiRep <- datMultiRep[datMultiRep$Multimodality_YN != "0", ]; datMultiRep <- datMultiRep[datMultiRep$Multimodality_YN != "NA", ]; datMultiRep <- subset(datMultiRep, !is.na(Multimodality_YN))
datMultiRep <- datMultiRep[!is.na(datMultiRep$SignalCombination), ]

# Omit all levels of Subject that contributed fewer than 5 multicomponent cases
nb.signals.per.ind <- data.frame(
  datMultiRep %>%
    group_by(Subject) %>%
    dplyr::summarize(count = n(), .groups = "drop")
)
colnames(nb.signals.per.ind) <- c("Subject", "NsignalsTot")
nb.signals.per.ind$Subject <- as.character(nb.signals.per.ind$Subject)

nb.signals.per.ind <- subset(nb.signals.per.ind, NsignalsTot > 4)
datMultiRep <- subset(datMultiRep, Subject %in% nb.signals.per.ind$Subject)

# Recompute after filtering individuals
nb.signals.per.ind <- data.frame(
  datMultiRep %>%
    group_by(Subject) %>%
    dplyr::summarize(count = n(), .groups = "drop")
)
colnames(nb.signals.per.ind) <- c("Subject", "NsignalsTot")
nb.signals.per.ind$Subject <- as.character(nb.signals.per.ind$Subject)

datMultiRep$Behavior <- datMultiRep$SignalCombination
Data <- datMultiRep

# ---------------------------------------------------------------------------------------------

# Individual repertoire similarity for groups:
library(future.apply)

# Use NULL for full data, or e.g. 50 for rarefaction to 50 signals per individual
NperInd <- 10 #min(nb.signals.per.ind$NsignalsTot)
Data$Subject <- as.character(Data$Subject)
Data$Behavior <- as.character(Data$Behavior)
Data$Group <- as.character(Data$Group)

plan(multisession, workers = 10)

nSubsample <- 100
nPerm <- 1000

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
resultsMR_G_WB <- future_lapply(X=1:nSubsample, future.seed=TRUE, FUN=function(s) {

  if(is.null(NperInd)) { Data.subsampled <- Data } else {
    valid_inds <- nb.signals.per.ind$Subject[nb.signals.per.ind$NsignalsTot >= NperInd]
    
    Data.tmp <- subset(Data, Subject %in% valid_inds)
    
    Data.subsampled <- Data.tmp %>%
      group_by(Subject, .drop = TRUE) %>%
      slice_sample(n = NperInd, replace = FALSE) %>%
      ungroup()
  }
  
  dat.ind.bhv2 <- data.frame(
    Data.subsampled %>%
      group_by(Subject, Behavior) %>%
      dplyr::summarize(count = n(), .groups = "drop")
  )
  colnames(dat.ind.bhv2) <- c("Subject","Behavior","NsignalsPbehavior")
  
  IndRepSize <- aggregate(data=Data.subsampled, Behavior ~ Subject, function(x) length(unique(x)))
  colnames(IndRepSize) <- c("Subject","RepertoireSize") 

  # Individual repertoire similarity (A-A/B-B, A-B): 
  subjects <- as.character(IndRepSize$Subject)
  
  bhv.list <- split(dat.ind.bhv2$Behavior, dat.ind.bhv2$Subject, drop=TRUE)
  bhv.list <- bhv.list[subjects]
  nsub <- length(subjects)
  
  group.vec <- Data.subsampled %>%
    distinct(Subject, Group) %>%
    tibble::deframe()
  
  group.labels <- group.vec[subjects]
  
  stopifnot(all(names(bhv.list) == subjects))
  stopifnot(all(names(group.labels) == subjects))
  
  # step 1 : calculate Dice coefficient based on formula: dc = (2 x number of behaviours two inds have in common)/(R of ind1 + R ind2) ---
  Dice <- matrix(NA, nsub, nsub, dimnames=list(subjects, subjects))
  
  for(i in 1:nsub){
    Ind1 <- unique(bhv.list[[i]])
    RInd1 <- length(Ind1)
    for(j in i:nsub){
      Ind2 <- unique(bhv.list[[j]])
      RInd2 <- length(Ind2)
      ovrlp <- length(intersect(Ind1, Ind2))
      Dice[i,j] <- (2 * ovrlp) / (RInd1 + RInd2)
      Dice[j,i] <- Dice[i,j]
    }
  }
  diag(Dice) <- NA
  lower_idx <- lower.tri(Dice)
  
  # step 2: define within- and between-group dyads ---
  wth.btw <- outer(
    group.labels,
    group.labels,
    FUN = function(x,y) { ifelse(x==y, "Within", "Between") }
  )
  diag(wth.btw) <- NA
  
  # Extract upper triangle only (unique dyads)
  upper_idx <- upper.tri(Dice, diag = FALSE)
  # Build dyad-level table
  WBdd <- data.frame(
    Ind1 = rownames(Dice)[row(Dice)[upper_idx]],
    Ind2 = colnames(Dice)[col(Dice)[upper_idx]],
    Dice = Dice[upper_idx],
    WithinBetween = wth.btw[upper_idx],
    stringsAsFactors = FALSE
  )
  # Remove NA values if needed
  WBdd <- WBdd[!is.na(WBdd$Dice), ]
  
  id_info <- Data.subsampled %>%
    distinct(Subject, Group)
  
  WBdd <- WBdd %>%
    left_join(id_info, by = c("Ind1" = "Subject")) %>%
    rename(Group.ind1 = Group) %>%
    left_join(id_info, by = c("Ind2" = "Subject")) %>%
    rename(Group.ind2 = Group)
  
  # step 3: get within and between group Dice coefficients ---
  Emp.MDiceWithin <- mean(Dice[lower_idx & wth.btw == "Within"], na.rm=TRUE)
  Emp.MDiceBetween <- mean(Dice[lower_idx & wth.btw == "Between"], na.rm=TRUE)
  EF.obs <- Emp.MDiceWithin - Emp.MDiceBetween
  
  # step 4: matrix permutation test ---
  dat.Perm <- numeric(nPerm)
  
  for(k in 1:nPerm) {
    perm.labels <- sample(group.labels, replace=FALSE)
    perm.wth.btw <- outer(
      perm.labels,
      perm.labels,
      FUN = function(x,y) ifelse(x==y, "Within", "Between")
    )
    diag(perm.wth.btw) <- NA
    
    MDiceWithin <- mean(Dice[lower_idx & perm.wth.btw=="Within"], na.rm=TRUE)
    MDiceBetween <- mean(Dice[lower_idx & perm.wth.btw=="Between"], na.rm=TRUE)
    
    dat.Perm[k] <- MDiceWithin - MDiceBetween
  }

  list(
    EF = EF.obs,
    pMRG = (sum(abs(dat.Perm) >= abs(EF.obs)) + 1) / (nPerm + 1),
    MW = Emp.MDiceWithin,
    MB = Emp.MDiceBetween,
    Dyads = WBdd
  )
})

EF_MRG <- sapply(resultsMR_G_WB, `[[`, "EF")
pMRG <- sapply(resultsMR_G_WB, `[[`, "pMRG")
MW_MRG <- sapply(resultsMR_G_WB, `[[`, "MW")
MB_MRG <- sapply(resultsMR_G_WB, `[[`, "MB")

all_dyadsMRG <- dplyr::bind_rows(
  lapply(seq_along(resultsMR_G_WB), function(i) {
    cbind(Subsample = i, resultsMR_G_WB[[i]]$Dyads)
  })
)
mean_dyadsMRG <- all_dyadsMRG %>%
  group_by(Ind1, Ind2, Group.ind1, Group.ind2, WithinBetween) %>%
  summarize(
    MeanDice = mean(Dice, na.rm = TRUE),
    SDDice = sd(Dice, na.rm = TRUE),
    Nd = n(),
    .groups = "drop"
  )

summary(EF_MRG)    			# Effect size
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.02260 0.03071 0.03450 0.03528 0.03786 0.05294 
sd(EF_MRG)					# Standard deviation of effect size across rarefactions
#0.006345141
summary(pMRG)    		# P-value
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000999 0.000999 0.000999 0.001459 0.000999 0.007992 
mean(pMRG < 0.05)    	# P-value < 0.05
#1
mean(MW_MRG)   				# Mean Dice within
#0.0796644
mean(MB_MRG)   				# Mean Dice between
#0.04438416
sd(MW_MRG)
#0.007349279
sd(MB_MRG)
#0.003935282


# ---------------------------------------------------------------------------------------------
# Repertoire similarity within groups (A-A, B-B, C2-C2, C4-C4, L-L, K-K) vs between groups ----
# ---------------------------------------------------------------------------------------------
# Individual repertoire similarity for groups:
library(future.apply)

groups <- c("A", "B", "C2", "C4", "L", "K")
# Use NULL for full data, or e.g. 10 for rarefaction to 10 signals per individual
# Or use min(nb.signals.per.ind$NsignalsTot) to select the minimum across individuals (i.e. 5)
NperInd <- 10 #min(nb.signals.per.ind$NsignalsTot)
Data$Subject <- as.character(Data$Subject)
Data$Behavior <- as.character(Data$Behavior)
Data$Group <- as.character(Data$Group)

plan(multisession, workers = 10)

nSubsample <- 100
nPerm <- 1000

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
resultsMR_G_W <- future_lapply(X=1:nSubsample, future.seed=TRUE, FUN=function(s) {

  if(is.null(NperInd)) { Data.subsampled <- Data } else {
    valid_inds <- nb.signals.per.ind$Subject[nb.signals.per.ind$NsignalsTot >= NperInd]
    
    Data.tmp <- subset(Data, Subject %in% valid_inds)
    
	# Only keep groups with at least two eligible individuals
    valid_groups <- Data.tmp %>%
      dplyr::distinct(Subject, Group) %>%
      dplyr::count(Group) %>%
      dplyr::filter(n >= 2) %>%
      dplyr::pull(Group)
	Data.tmp <- subset(Data.tmp, Group %in% valid_groups)
	
    Data.subsampled <- Data.tmp %>%
      group_by(Subject, .drop = TRUE) %>%
      slice_sample(n = NperInd, replace = FALSE) %>%
      ungroup()
  }
  
  dat.ind.bhv2 <- data.frame(
    Data.subsampled %>%
      group_by(Subject, Behavior) %>%
      dplyr::summarize(count = n(), .groups = "drop")
  )
  colnames(dat.ind.bhv2) <- c("Subject","Behavior","NsignalsPbehavior")
  
  IndRepSize <- aggregate(data=Data.subsampled, Behavior ~ Subject, function(x) length(unique(x)))
  colnames(IndRepSize) <- c("Subject","RepertoireSize") 

  # Individual repertoire similarity (A-A/B-B, A-B): 
  subjects <- as.character(IndRepSize$Subject)
  
  bhv.list <- split(dat.ind.bhv2$Behavior, dat.ind.bhv2$Subject, drop=TRUE)
  bhv.list <- bhv.list[subjects]
  nsub <- length(subjects)
  
  group.vec <- Data.subsampled %>%
    distinct(Subject, Group) %>%
    tibble::deframe()
  
  group.labels <- group.vec[subjects]
  groups.use <- groups[groups %in% unique(group.labels)]
  
  stopifnot(all(names(bhv.list) == subjects))
  stopifnot(all(names(group.labels) == subjects))
  
  # step 1 : calculate Dice coefficient based on formula: dc = (2 x number of behaviours two inds have in common)/(R of ind1 + R ind2) ---
  Dice <- matrix(NA, nsub, nsub, dimnames=list(subjects, subjects))
  
  for(i in 1:nsub){
    Ind1 <- unique(bhv.list[[i]])
    RInd1 <- length(Ind1)
    for(j in i:nsub){
      Ind2 <- unique(bhv.list[[j]])
      RInd2 <- length(Ind2)
      ovrlp <- length(intersect(Ind1, Ind2))
      Dice[i,j] <- (2 * ovrlp) / (RInd1 + RInd2)
      Dice[j,i] <- Dice[i,j]
    }
  }
  diag(Dice) <- NA
  lower_idx <- lower.tri(Dice)

  # Matrix saying which within-group dyad each pair belongs to
  wth.btw2 <- outer(
    group.labels,
    group.labels,
    FUN = function(x, y) ifelse(x == y, x, NA_character_))
  diag(wth.btw2) <- NA_character_
  
  # Empirical group-specific mean Dice values
  mu.emp <- sapply(groups.use, function(g) {
    mean(Dice[lower_idx & wth.btw2 == g], na.rm = TRUE)
  })

  names(mu.emp) <- groups.use
  stopifnot(!any(is.na(mu.emp)))
  
  # Weights = number of within-group dyads per group
  group_sizes <- table(factor(group.labels, levels = groups.use))
  w <- group_sizes * (group_sizes - 1) / 2
  stopifnot(all(groups.use %in% unique(group.labels)))
  stopifnot(all(w > 0))
  
  # Weighted omnibus statistic
  mean.w.emp <- weighted.mean(mu.emp, w = w, na.rm = TRUE)
  SS.emp <- sum(w * (mu.emp - mean.w.emp)^2, na.rm = TRUE)
  
  # Pairwise empirical differences
  contrast.mat <- combn(groups.use, 2)

  pairwise.emp <- apply(contrast.mat, 2, function(z) {
    mu.emp[z[1]] - mu.emp[z[2]]
  })

  contrast.names <- apply(contrast.mat, 2, paste, collapse = "-")
  names(pairwise.emp) <- contrast.names
  
  # Permutation test
  dat.Perm2 <- numeric(nPerm)
  perm.pairwise <- matrix(NA_real_, nrow = nPerm, ncol = length(pairwise.emp))
  colnames(perm.pairwise) <- contrast.names
  
  for(k in 1:nPerm) {
    
    perm.labels <- sample(group.labels, replace = FALSE)
    
    perm.wth.btw2 <- outer(
      perm.labels,
      perm.labels,
      FUN = function(x, y) ifelse(x == y, x, NA_character_))
    diag(perm.wth.btw2) <- NA_character_
    
    mu <- sapply(groups.use, function(g) {
      mean(Dice[lower_idx & perm.wth.btw2 == g], na.rm = TRUE)
    })
    names(mu) <- groups.use
    stopifnot(!any(is.na(mu)))

	mean.w <- weighted.mean(mu, w = w, na.rm = TRUE)
    dat.Perm2[k] <- sum(w * (mu - mean.w)^2, na.rm = TRUE)
    
	perm.pairwise[k, ] <- apply(contrast.mat, 2, function(z) {
      mu[z[1]] - mu[z[2]]
    })
  }
  
  p.omnibus <- (sum(dat.Perm2 >= SS.emp) + 1) / (nPerm + 1)
  
  p.pairwise <- sapply(seq_along(pairwise.emp), function(i) {
    valid <- !is.na(perm.pairwise[, i])
    (sum(abs(perm.pairwise[valid, i]) >= abs(pairwise.emp[i])) + 1) /
      (sum(valid) + 1)
  })
  names(p.pairwise) <- contrast.names
  
  p.pairwise.holm <- p.adjust(p.pairwise, method = "holm")
  names(p.pairwise.holm) <- contrast.names
  
  list(
    SS = SS.emp,
    p.omnibus = p.omnibus,
    mu.emp = mu.emp,
    pairwise.emp = pairwise.emp,
    p.pairwise = p.pairwise,
    p.pairwise.holm = p.pairwise.holm
  )
})

SS_MRG <- sapply(resultsMR_G_W, `[[`, "SS")
p.omnibus_MRG <- sapply(resultsMR_G_W, `[[`, "p.omnibus")
mu.emp_MRG <- sapply(resultsMR_G_W, `[[`, "mu.emp")
pairwise.emp_MRG <- sapply(resultsMR_G_W, `[[`, "pairwise.emp")
p.pairwise_MRG <- sapply(resultsMR_G_W, `[[`, "p.pairwise")
p.pairwise.holm_MRG <- sapply(resultsMR_G_W, `[[`, "p.pairwise.holm")

summary(SS_MRG)    						# How stable are the omnibus effect size is across subsampling?
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.02920 0.06455 0.08214 0.08739 0.10860 0.19344 
sd(SS_MRG)								# Standard deviation of the Sum of Squares
#0.03481935
summary(p.omnibus_MRG)    				# omnibus P-values
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.005994 0.090659 0.139361 0.178541 0.227522 0.586414 
mean(p.omnibus_MRG < 0.05)  			# Does the omnibus significance remain robust across subsampling?
#0.1
rowMeans(mu.emp_MRG)    				# Average repertoire similarity and its variability per group
#         A         C2          L          K 
#0.09586526 0.08355756 0.03596343 0.06133126 
apply(mu.emp_MRG, 1, sd)
#          A          C2           L           K 
#0.017018480 0.010278513 0.005830747 0.017683520 
rowMeans(pairwise.emp_MRG)
#      A-C2         A-L         A-K        C2-L        C2-K         L-K 
#0.01230769  0.05990182  0.03453399  0.04759413  0.02222630 -0.02536783 
apply(pairwise.emp_MRG, 1, sd)
#      A-C2        A-L        A-K       C2-L       C2-K        L-K 
#0.01953346 0.01819098 0.02474750 0.01150121 0.01978767 0.01912151 
rowMeans(p.pairwise_MRG)
#      A-C2        A-L        A-K       C2-L       C2-K        L-K 
#0.49822178 0.09955045 0.41585415 0.08946054 0.50660340 0.52708292 
rowMeans(p.pairwise_MRG < 0.05)			# In what fraction of rarefaction replicates was this contrast significant?
#A-C2  A-L  A-K C2-L C2-K  L-K 
#0.03 0.49 0.05 0.34 0.00 0.00 
rowMeans(p.pairwise.holm_MRG < 0.05)	# Which contrasts remain significant after multiple-testing correction consistently across rarefactions?
#A-C2  A-L  A-K C2-L C2-K  L-K 
#0.02 0.04 0.02 0.00 0.00 0.00 


# ---------------------------------------------------------------------------------------------
# Boxplot MeanDice within/between (steps a) --------------------------------------------------
# ---------------------------------------------------------------------------------------------
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

levels(mean_dyadsMRG$Group.ind1) <- c("A", "B", "L", "C2", "K") # dice = ddMulti ?

WBddMultiWithin <- subset(mean_dyadsMRG, mean_dyadsMRG$WithinBetween=="Within")

F3 <- ggplot() + 
  geom_boxplot(WBddMultiWithin, mapping=aes(x = Group.ind1, y = MeanDice), width = 0.9, fill = c("#007ABB","#259C39", "#E7B800", "#ba0f09")) +
  geom_boxplot(mean_dyadsMRG, mapping=aes(x = WithinBetween, y = MeanDice), width = 0.9, fill = c("grey90","grey90")) +
  geom_point(WBddMultiWithin, mapping=aes(x = Group.ind1, y = MeanDice), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_point(mean_dyadsMRG, mapping=aes(x = WithinBetween, y = MeanDice), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey40") + 
  theme_angele_ss +
  ggtitle("(a)") +
  scale_y_continuous("Repertoire similarity among individuals") +
  scale_x_discrete(" ",
                   limits = c("A", "L", "C2", "K", "Between", "Within"),
                   labels = c("Within \ngroup A", "Within \ngroup L", "Within \ngroup C2", "Within \ngroup K", "Between \ngroups", "Within \ngroups"))+
  #facet_wrap(~MeanDice, scales='free_x')+
  stat_summary(WBddMultiWithin, mapping=aes(x = Group.ind1, y = MeanDice), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3) +
  stat_summary(mean_dyadsMRG, mapping=aes(x = WithinBetween, y = MeanDice), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3)

# + theme(axis.text.x = element_blank()) 

print(F3)
