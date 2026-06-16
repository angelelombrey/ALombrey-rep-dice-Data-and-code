## 4. GROUP SIGNAL REPERTOIRE SIMILARITY ##

# PACKAGES
library(lme4);library(tidyverse); library(dplyr);library(writexl);library(ggpubr) 
library(tidyr);library(car);library(RColorBrewer);library(ggplot2) #library(plyr);
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
                         title = element_text(size = 18, family="sans"),
                         strip.text = element_text(size = 15))
# PREPARE DATA
path <- "D:\\Documents\\CoAuthoredManuscripts\\2026__Lombrey__Chimpanzee_Plasticity\\DataAnalysis\\"
library(openxlsx)
Data <- read.xlsx(xlsxFile=paste0(path,"Data.xlsx"), sheet=1, colNames=TRUE)

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
nb.signals.per.ind <- data.frame(Data %>% group_by(Subject) %>% dplyr::summarize(count=n(), .groups = "drop")) ; colnames(nb.signals.per.ind) <- c("Subject","NsignalsTot") 
nb.signals.per.ind$Subject <- as.character(nb.signals.per.ind$Subject)
nb.signals.per.ind <- subset(nb.signals.per.ind, nb.signals.per.ind$NsignalsTot>=30) # Keeps only the values in "column" that appear more than 30 times (N >= 30)
Data <- subset(Data, Subject %in% nb.signals.per.ind$Subject)

Data <- Data[!is.na(Data$Group), ]
nb.signals.per.ind <- data.frame(Data %>% group_by(Subject) %>% dplyr::summarize(count=n(), .groups = "drop")) ; colnames(nb.signals.per.ind) <- c("Subject","NsignalsTot") 
nb.signals.per.ind$Subject <- as.character(nb.signals.per.ind$Subject)
nb.signals.per.ind <- subset(nb.signals.per.ind, NsignalsTot >= 30)
Data <- subset(Data, Subject %in% nb.signals.per.ind$Subject)

# Sampling effort
Datareduit <- Data[Data$Duplicated != "T", ]  
SamplingEffort <- data.frame(Datareduit %>% group_by(Subject) %>% dplyr::summarize(count=n())) ; colnames(SamplingEffort) <- c("Subject","TotComActs") 
Data <- merge(SamplingEffort, Data, by=c("Subject"))

IndRepSize <- aggregate(data=Data, Behavior ~ Subject, function(x) length(unique(x))) ; colnames(IndRepSize) <- c("Subject","RepertoireSize") 
dat.ind.bhv2 <- data.frame(Data %>% group_by(Subject,Behavior) %>% dplyr::summarize(count=n(), .groups = "drop")) ; colnames(dat.ind.bhv2) <- c("Subject","Behavior","NsignalsPbehavior") 

# ---------------------------------------------------------------------------------------------

# Individual repertoire similarity for groups:
library(future.apply)

# Use NULL for full data, or e.g. 50 for rarefaction to 50 signals per individual
NperInd <- 50 # OR: min(nb.signals.per.ind$NsignalsTot)
Data$Subject <- as.character(Data$Subject)
Data$Behavior <- as.character(Data$Behavior)
Data$Group <- as.character(Data$Group)

plan(multisession, workers = 10)

nSubsample <- 100
nPerm <- 1000

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
resultsSR_G_WB <- future_lapply(X=1:nSubsample, future.seed=TRUE, FUN=function(s) {

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
    pSRG = (sum(abs(dat.Perm) >= abs(EF.obs)) + 1) / (nPerm + 1),
    MW = Emp.MDiceWithin,
    MB = Emp.MDiceBetween,
    Dyads = WBdd
  )
})

EF_SRG <- sapply(resultsSR_G_WB, `[[`, "EF")
pSRG <- sapply(resultsSR_G_WB, `[[`, "pSRG")
MW_SRG <- sapply(resultsSR_G_WB, `[[`, "MW")
MB_SRG <- sapply(resultsSR_G_WB, `[[`, "MB")

all_dyadsSRG <- dplyr::bind_rows(
  lapply(seq_along(resultsSR_G_WB), function(i) {
    cbind(Subsample = i, resultsSR_G_WB[[i]]$Dyads)
  })
)
mean_dyadsSRG <- all_dyadsSRG %>%
  group_by(Ind1, Ind2, Group.ind1, Group.ind2, WithinBetween) %>%
  summarize(
    MeanDice = mean(Dice, na.rm = TRUE),
    SDDice = sd(Dice, na.rm = TRUE),
    Nd = n(),
    .groups = "drop"
  )

summary(EF_SRG)    			# Effect size
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.07426 0.08372 0.08826 0.08773 0.09129 0.10356 
sd(EF_SRG)					# Standard deviation of effect size across rarefactions
#0.005800659
summary(pSRG)    		# P-value
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000999 0.000999 0.000999 0.000999 0.000999 0.000999 
mean(pSRG < 0.05)    	# P-value < 0.05
mean(MW_SRG)   				# Mean Dice within
#0.5649512
mean(MB_SRG)   				# Mean Dice between
#0.4772214
sd(MW_SRG)
#0.006036146
sd(MB_SRG)
#0.004795485

# ---------------------------------------------------------------------------------------------
# Repertoire similarity within groups (A-A, B-B, C2-C2, C4-C4, L-L, K-K) vs between groups ----
# ---------------------------------------------------------------------------------------------

groups <- c("A", "B", "C2", "C4", "L", "K")
# Use NULL for full data, or e.g. 50 for rarefaction to 50 signals per individual
# Or use min(nb.signals.per.ind$NsignalsTot) to select the minimum across individuals (i.e. 30)
NperInd <- 50 # min(nb.signals.per.ind$NsignalsTot)
Data$Subject <- as.character(Data$Subject)
Data$Behavior <- as.character(Data$Behavior)
Data$Group <- as.character(Data$Group)

plan(multisession, workers = 10)

nSubsample <- 100
nPerm <- 1000

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
resultsSR_G_W <- future_lapply(X=1:nSubsample, future.seed=TRUE, FUN=function(s) {

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

  # Matrix saying which within-group dyad each pair belongs to
  wth.btw2 <- outer(
    group.labels,
    group.labels,
    FUN = function(x, y) ifelse(x == y, x, NA_character_))
  diag(wth.btw2) <- NA_character_
  
  # Empirical group-specific mean Dice values
  mu.emp <- sapply(groups, function(g) {
    mean(Dice[lower_idx & wth.btw2 == g], na.rm = TRUE)
  })

  names(mu.emp) <- groups
  stopifnot(!any(is.na(mu.emp)))
  
  # Weights = number of within-group dyads per group
  group_sizes <- table(factor(group.labels, levels = groups))
  w <- group_sizes * (group_sizes - 1) / 2
  stopifnot(all(groups %in% unique(group.labels)))
  stopifnot(all(w > 0))
  
  # Weighted omnibus statistic
  mean.w.emp <- weighted.mean(mu.emp, w = w, na.rm = TRUE)
  SS.emp <- sum(w * (mu.emp - mean.w.emp)^2, na.rm = TRUE)
  
  # Pairwise empirical differences
  pairwise.emp <- c(
    mu.emp["A"]  - mu.emp["B"],
    mu.emp["A"]  - mu.emp["C2"],
    mu.emp["A"]  - mu.emp["C4"],
    mu.emp["A"]  - mu.emp["L"],
    mu.emp["A"]  - mu.emp["K"],
    mu.emp["B"]  - mu.emp["C2"],
    mu.emp["B"]  - mu.emp["C4"],
    mu.emp["B"]  - mu.emp["L"],
    mu.emp["B"]  - mu.emp["K"],
    mu.emp["C2"] - mu.emp["C4"],
    mu.emp["C2"] - mu.emp["L"],
    mu.emp["C2"] - mu.emp["K"],
    mu.emp["C4"] - mu.emp["L"],
    mu.emp["C4"] - mu.emp["K"],
    mu.emp["L"]  - mu.emp["K"]
  )
  
  contrast.names <- c(
    "A-B","A-C2","A-C4","A-L","A-K",
    "B-C2","B-C4","B-L","B-K",
    "C2-C4","C2-L","C2-K",
    "C4-L","C4-K",
	"L-K"
  )
  
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
    
    mu <- sapply(groups, function(g) {
      mean(Dice[lower_idx & perm.wth.btw2 == g], na.rm = TRUE)
    })
    names(mu) <- groups
    stopifnot(!any(is.na(mu)))

	mean.w <- weighted.mean(mu, w = w, na.rm = TRUE)
    dat.Perm2[k] <- sum(w * (mu - mean.w)^2, na.rm = TRUE)
    
    perm.pairwise[k, ] <- c(
      mu["A"]  - mu["B"],
      mu["A"]  - mu["C2"],
      mu["A"]  - mu["C4"],
      mu["A"]  - mu["L"],
      mu["A"]  - mu["K"],
      mu["B"]  - mu["C2"],
      mu["B"]  - mu["C4"],
      mu["B"]  - mu["L"],
      mu["B"]  - mu["K"],
      mu["C2"] - mu["C4"],
      mu["C2"] - mu["L"],
      mu["C2"] - mu["K"],
      mu["C4"] - mu["L"],
      mu["C4"] - mu["K"],
      mu["L"]  - mu["K"]
    )
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

SS_SRG <- sapply(resultsSR_G_W, `[[`, "SS")
p.omnibus_SRG <- sapply(resultsSR_G_W, `[[`, "p.omnibus")
mu.emp_SRG <- sapply(resultsSR_G_W, `[[`, "mu.emp")
pairwise.emp_SRG <- sapply(resultsSR_G_W, `[[`, "pairwise.emp")
p.pairwise_SRG <- sapply(resultsSR_G_W, `[[`, "p.pairwise")
p.pairwise.holm_SRG <- sapply(resultsSR_G_W, `[[`, "p.pairwise.holm")

summary(SS_SRG)    						# How stable are the omnibus effect size is across subsampling?
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3637  0.6315  0.7336  0.7481  0.8379  1.5208 
sd(SS_SRG)								# Standard deviation of the Sum of Squares
#0.174949
summary(p.omnibus_SRG)    				# omnibus P-values
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.003996 0.066933 0.106394 0.113966 0.147852 0.440559 
mean(p.omnibus_SRG < 0.05)  			# Does the omnibus significance remain robust across subsampling?
#0.14
rowMeans(mu.emp_SRG)    				# Average repertoire similarity and its variability per group
#        A         B        C2        C4         L         K 
#0.5789804 0.2869004 0.5755533 0.6348287 0.4903216 0.5878722 
apply(mu.emp_SRG, 1, sd)
#          A           B          C2          C4           L           K 
#0.013696109 0.045055725 0.009247084 0.030948479 0.008961824 0.010816373 
rowMeans(pairwise.emp_SRG)
#        A-B         A-C2         A-C4          A-L          A-K         B-C2         B-C4          B-L          B-K        C2-C4         C2-L         C2-K         C4-L         C4-K          L-K 
#0.292079965  0.003427063 -0.055848371  0.088658775 -0.008891803 -0.288652902 -0.347928337 -0.203421190 -0.300971768 -0.059275435  0.085231712 -0.012318866  0.144507147  0.046956569 -0.097550578 
apply(pairwise.emp_SRG, 1, sd)
#       A-B       A-C2       A-C4        A-L        A-K       B-C2       B-C4        B-L        B-K      C2-C4       C2-L       C2-K       C4-L       C4-K        L-K 
#0.04848537 0.01747612 0.03405675 0.01723218 0.01873755 0.04595276 0.05340338 0.04391660 0.04525622 0.03148145 0.01246937 0.01411194 0.03455420 0.03380299 0.01444464 
rowMeans(p.pairwise_SRG)
#       A-B       A-C2       A-C4        A-L        A-K       B-C2       B-C4        B-L        B-K      C2-C4       C2-L       C2-K       C4-L       C4-K        L-K 
#0.03800200 0.76352647 0.52624376 0.10199800 0.75174825 0.03453546 0.02559441 0.14700300 0.02994006 0.47205794 0.04964036 0.72031968 0.09764236 0.57566434 0.05759241 
rowMeans(p.pairwise_SRG < 0.05)			# In what fraction of rarefaction replicates was this contrast significant?
# A-B  A-C2  A-C4   A-L   A-K  B-C2  B-C4   B-L   B-K C2-C4  C2-L  C2-K  C4-L  C4-K   L-K 
#0.71  0.00  0.00  0.19  0.00  0.69  0.85  0.08  0.86  0.00  0.59  0.00  0.37  0.00  0.51 
rowMeans(p.pairwise.holm_SRG < 0.05)	# Which contrasts remain significant after multiple-testing correction consistently across rarefactions?
# A-B  A-C2  A-C4   A-L   A-K  B-C2  B-C4   B-L   B-K C2-C4  C2-L  C2-K  C4-L  C4-K   L-K 
#0.05  0.00  0.00  0.00  0.00  0.05  0.08  0.00  0.04  0.00  0.02  0.00  0.01  0.00  0.00 


# ---------------------------------------------------------------------------------------------
# Boxplot Dice within/between (steps a) -------------------------------------------------------
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
                         axis.title.y = element_text(size = 19, vjust = 2, family="sans"),
                         legend.text=  element_text(size = 15, family="sans", margin = margin(t = 10)),
                         legend.title = element_text(size = 15, vjust = 2, family="sans"),
                         legend.key = element_blank(),
                         legend.position = "right",
                         legend.spacing.x = unit(0.2, 'cm'),
                         title = element_text(size = 20, family="sans"),
                         strip.text = element_text(size = 15))

levels(mean_dyadsSRG$Group.ind1) <- c("A", "B", "L", "C2", "C4", "K") # dice = dd ?

WBddWithin <- subset(mean_dyadsSRG, mean_dyadsSRG$WithinBetween=="Within")


F1 <- ggplot() + 
  geom_boxplot(WBddWithin, mapping=aes(x = Group.ind1, y = MeanDice), width = 0.9, fill = c("#007ABB", "#00AFBB","#259C39", "#E7B800","#E79A00", "#ba0f09")) +
  geom_boxplot(mean_dyadsSRG, mapping=aes(x = WithinBetween, y = MeanDice), width = 0.9, fill = c("grey90","grey90")) +
  geom_point(WBddWithin, mapping=aes(x = Group.ind1, y = MeanDice), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_point(mean_dyadsSRG, mapping=aes(x = WithinBetween, y = MeanDice), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "grey40") + 
  theme_angele_ss +
  ggtitle("(a)") +
  scale_y_continuous("Repertoire similarity among individuals", breaks = seq(0, 1, by = 0.2)) +
  scale_x_discrete(" ",
                   limits = c("A", "B", "L", "C2", "C4", "K", "Between", "Within"),
                   labels = c("Within \ngroup A", "Within \ngroup B", "Within \ngroup L", "Within \ngroup C2", "Within \ngroup C4", "Within \ngroup K", "Between \ngroups", "Within \ngroups"))+
  #facet_wrap(~Dice, scales='free_x')+
  stat_summary(WBddWithin, mapping=aes(x = Group.ind1, y = MeanDice), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3) +
  stat_summary(mean_dyadsSRG, mapping=aes(x = WithinBetween, y = MeanDice), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3)

# + theme(axis.text.x = element_blank()) 

print(F1)
