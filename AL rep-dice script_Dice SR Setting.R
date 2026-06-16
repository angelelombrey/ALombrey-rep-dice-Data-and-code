## 5. SETTING SIGNAL REPERTOIRE SIMILARITY ##

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

Data <- Data[!is.na(Data$Setting), ]
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

# Individual repertoire similarity for settings (Captive-Captive/Sanctuary-Sanctuary, Captive-Sanctuary): 
library(future.apply)

# Use NULL for full data, or e.g. 50 for rarefaction to 50 signals per individual
NperInd <- 50 #min(nb.signals.per.ind$NsignalsTot)
Data$Subject <- as.character(Data$Subject)
Data$Behavior <- as.character(Data$Behavior)
Data$Setting <- as.character(Data$Setting)

plan(multisession, workers = 10)

nSubsample <- 100
nPerm <- 1000

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
resultsSR_S_WB <- future_lapply(X=1:nSubsample, future.seed=TRUE, FUN=function(s) {

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

  # Individual repertoire similarity among settings: 
  subjects <- as.character(IndRepSize$Subject)
  
  bhv.list <- split(dat.ind.bhv2$Behavior, dat.ind.bhv2$Subject, drop=TRUE)
  bhv.list <- bhv.list[subjects]
  nsub <- length(subjects)
  
  setting.vec <- Data.subsampled %>%
    distinct(Subject, Setting) %>%
    tibble::deframe()
  
  setting.labels <- setting.vec[subjects]
  
  stopifnot(all(names(bhv.list) == subjects))
  stopifnot(all(names(setting.labels) == subjects))
  
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
  
  # step 2: define within- and between-setting dyads ---
  wth.btw <- outer(
    setting.labels,
    setting.labels,
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
    distinct(Subject, Setting)
  
  WBdd <- WBdd %>%
    left_join(id_info, by = c("Ind1" = "Subject")) %>%
    rename(Setting.ind1 = Setting) %>%
    left_join(id_info, by = c("Ind2" = "Subject")) %>%
    rename(Setting.ind2 = Setting)
  
  # step 3: get within and between setting Dice coefficients ---
  Emp.MDiceWithin <- mean(Dice[lower_idx & wth.btw == "Within"], na.rm=TRUE)
  Emp.MDiceBetween <- mean(Dice[lower_idx & wth.btw == "Between"], na.rm=TRUE)
  EF.obs <- Emp.MDiceWithin - Emp.MDiceBetween
  
  # step 4: matrix permutation test ---
  dat.Perm <- numeric(nPerm)
  
  for(k in 1:nPerm) {
    perm.labels <- sample(setting.labels, replace=FALSE)
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
    pSRS = (sum(abs(dat.Perm) >= abs(EF.obs)) + 1) / (nPerm + 1),
    MW = Emp.MDiceWithin,
    MB = Emp.MDiceBetween,
    Dyads = WBdd
  )
})

EF_SRS <- sapply(resultsSR_S_WB, `[[`, "EF")
pSRS <- sapply(resultsSR_S_WB, `[[`, "pSRS")
MW_SRS <- sapply(resultsSR_S_WB, `[[`, "MW")
MB_SRS <- sapply(resultsSR_S_WB, `[[`, "MB")

all_dyadsSRS <- dplyr::bind_rows(
  lapply(seq_along(resultsSR_S_WB), function(i) {
    cbind(Subsample = i, resultsSR_S_WB[[i]]$Dyads)
  })
)
mean_dyadsSRS <- all_dyadsSRS %>%
  group_by(Ind1, Ind2, Setting.ind1, Setting.ind2, WithinBetween) %>%
  summarize(
    MeanDice = mean(Dice, na.rm = TRUE),
    SDDice = sd(Dice, na.rm = TRUE),
    Nd = n(),
    .groups = "drop"
  )

summary(EF_SRS)    			# Effect size
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.04281 0.05379 0.05672 0.05657 0.05877 0.07015 
sd(EF_SRS)					# Standard deviation of effect size across rarefactions
#0.004242743
summary(pSRS)    		# P-value
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000999 0.000999 0.000999 0.000999 0.000999 0.000999 
mean(pSRS < 0.05)    	# P-value < 0.05
#1
mean(MW_SRS)   				# Mean Dice within
#0.5334221
mean(MB_SRS)   				# Mean Dice between
#0.4768491
sd(MW_SRS)
#0.005446644
sd(MB_SRS)
#0.004862747


# ---------------------------------------------------------------------------------------------
# Repertoire similarity within setting Captive (C-C), within setting Sanctuary (S-S) , within setting Wild (W-W)
# ---------------------------------------------------------------------------------------------
settings <- c("Captive", "Sanctuary", "Wild")
# Use NULL for full data, or e.g. 50 for rarefaction to 50 signals per individual
# Or use min(nb.signals.per.ind$NsignalsTot) to select the minimum across individuals (i.e. 30)
NperInd <- 50 #min(nb.signals.per.ind$NsignalsTot)
Data$Subject <- as.character(Data$Subject)
Data$Behavior <- as.character(Data$Behavior)
Data$Setting <- as.character(Data$Setting)

plan(multisession, workers = 10)

nSubsample <- 100
nPerm <- 1000

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
resultsSR_S_WB <- future_lapply(X=1:nSubsample, future.seed=TRUE, FUN=function(s) {

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
  
  setting.vec <- Data.subsampled %>%
    distinct(Subject, Setting) %>%
    tibble::deframe()
  
  setting.labels <- setting.vec[subjects]
  
  stopifnot(all(names(bhv.list) == subjects))
  stopifnot(all(names(setting.labels) == subjects))
  
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

  # Matrix saying which within-setting dyad each pair belongs to
  wth.btw2 <- outer(
    setting.labels,
    setting.labels,
    FUN = function(x, y) ifelse(x == y, x, NA_character_))
  diag(wth.btw2) <- NA_character_
  
  # Empirical setting-specific mean Dice values
  mu.emp <- sapply(settings, function(g) {
    mean(Dice[lower_idx & wth.btw2 == g], na.rm = TRUE)
  })

  names(mu.emp) <- settings
  stopifnot(!any(is.na(mu.emp)))
  
  # Weights = number of within-setting dyads per setting
  setting_sizes <- table(factor(setting.labels, levels = settings))
  w <- setting_sizes * (setting_sizes - 1) / 2
  stopifnot(all(settings %in% unique(setting.labels)))
  stopifnot(all(w > 0))
  
  # Weighted omnibus statistic
  mean.w.emp <- weighted.mean(mu.emp, w = w, na.rm = TRUE)
  SS.emp <- sum(w * (mu.emp - mean.w.emp)^2, na.rm = TRUE)
  
  # Pairwise empirical differences
  pairwise.emp <- c(
    mu.emp["Captive"]  - mu.emp["Sanctuary"],
    mu.emp["Captive"]  - mu.emp["Wild"],
    mu.emp["Sanctuary"]  - mu.emp["Wild"]
  )
  
  contrast.names <- c("Captive-Sanctuary","Captive-Wild","Sanctuary-Wild")
  
  names(pairwise.emp) <- contrast.names
  
  # Permutation test
  dat.Perm2 <- numeric(nPerm)
  perm.pairwise <- matrix(NA_real_, nrow = nPerm, ncol = length(pairwise.emp))
  colnames(perm.pairwise) <- contrast.names
  
  for(k in 1:nPerm) {
    
    perm.labels <- sample(setting.labels, replace = FALSE)
    
    perm.wth.btw2 <- outer(
      perm.labels,
      perm.labels,
      FUN = function(x, y) ifelse(x == y, x, NA_character_))
    diag(perm.wth.btw2) <- NA_character_
    
    mu <- sapply(settings, function(g) {
      mean(Dice[lower_idx & perm.wth.btw2 == g], na.rm = TRUE)
    })
    names(mu) <- settings
    stopifnot(!any(is.na(mu)))

	mean.w <- weighted.mean(mu, w = w, na.rm = TRUE)
    dat.Perm2[k] <- sum(w * (mu - mean.w)^2, na.rm = TRUE)
    
    perm.pairwise[k, ] <- c(
      mu["Captive"]  - mu["Sanctuary"],
      mu["Captive"]  - mu["Wild"],
      mu["Sanctuary"]  - mu["Wild"]
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

SS_SRS <- sapply(resultsSR_S_WB, `[[`, "SS")
p.omnibus_SRS <- sapply(resultsSR_S_WB, `[[`, "p.omnibus")
mu.emp_SRS <- sapply(resultsSR_S_WB, `[[`, "mu.emp")
pairwise.emp_SRS <- sapply(resultsSR_S_WB, `[[`, "pairwise.emp")
p.pairwise_SRS <- sapply(resultsSR_S_WB, `[[`, "p.pairwise")
p.pairwise.holm_SRS <- sapply(resultsSR_S_WB, `[[`, "p.pairwise.holm")

summary(SS_SRS)    						# How stable are the omnibus effect size is across subsampling?
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.2220  0.5389  0.6670  0.6906  0.8295  1.2137 
sd(SS_SRS)								# Standard deviation of the Sum of Squares
#0.2145786
summary(p.omnibus_SRS)    				# omnibus P-values
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.02298 0.07992 0.12787 0.14789 0.19880 0.46454 
mean(p.omnibus_SRS < 0.05)  			# Does the omnibus significance remain robust across subsampling?
#0.09
rowMeans(mu.emp_SRS)    				# Average repertoire similarity and its variability per group
#  Captive Sanctuary      Wild 
#0.5029355 0.5482853 0.5872770 
apply(mu.emp_SRS, 1, sd)
#    Captive   Sanctuary        Wild 
#0.006509037 0.008486232 0.011224118 
rowMeans(pairwise.emp_SRS)
#Captive-Sanctuary      Captive-Wild    Sanctuary-Wild 
#      -0.04534982       -0.08434146       -0.03899164 
apply(pairwise.emp_SRS, 1, sd)
#Captive-Sanctuary      Captive-Wild    Sanctuary-Wild 
#       0.01127220        0.01304797        0.01361987 
rowMeans(p.pairwise_SRS)
#Captive-Sanctuary      Captive-Wild    Sanctuary-Wild 
#       0.19598402        0.06178821        0.38064935 
rowMeans(p.pairwise_SRS < 0.05)			# In what fraction of rarefaction replicates was this contrast significant?
#Captive-Sanctuary      Captive-Wild    Sanctuary-Wild 
#             0.04              0.54              0.00 
rowMeans(p.pairwise.holm_SRS < 0.05)	# Which contrasts remain significant after multiple-testing correction consistently across rarefactions?
#Captive-Sanctuary      Captive-Wild    Sanctuary-Wild 
#              0.0               0.1               0.0 


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

levels(mean_dyadsSRS$Setting.ind1) <- c("Captive", "Sanctuary", "Wild") # Dice = dd ?

WBddSETWithin <- subset(mean_dyadsSRS, mean_dyadsSRS$WithinBetween=="Within")

F2 <- ggplot() + 
  geom_boxplot(WBddSETWithin, mapping=aes(x = Setting.ind1, y = MeanDice), width = 0.9, fill = c("#068591", "#ba8211", "#ba0f09")) +
  geom_boxplot(mean_dyadsSRS, mapping=aes(x = WithinBetween, y = MeanDice), width = 0.9, fill = c("grey90","grey90")) +
  geom_point(WBddSETWithin, mapping=aes(x = Setting.ind1, y = MeanDice), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_point(mean_dyadsSRS, mapping=aes(x = WithinBetween, y = MeanDice), position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey40") + 
  theme_angele_ss +
  ggtitle("(b)") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_x_discrete(" ",
                   limits = c("Captive", "Sanctuary", "Wild", "Between", "Within"),
                   labels = c("Within \nZoo", "Within \nSanctuary", "Within \nWild", "Between \nsettings", "Within \nsettings"))+
  #facet_wrap(~Dice, scales='free_x')+
  stat_summary(WBddSETWithin, mapping=aes(x = Setting.ind1, y = MeanDice), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3) +
  stat_summary(mean_dyadsSRS, mapping=aes(x = WithinBetween, y = MeanDice), fun=mean, geom="point",shape =23, fill ="black",position=position_dodge(.9), 
               color="black", size=3) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

print(F2)
ggarrange(F1, F2, widths = c(1.7,1))

# Setting VS Group comparison
Dice_group_set <- data.frame(
  value = c(MW_SRG, MW_SRS),
  group = rep(c("MW_SRG", "MW_SRS"),
              c(length(MW_SRG), length(MW_SRS)))
)
DiceB_group_set <- data.frame(
  value = c(MB_SRG, MB_SRS),
  group = rep(c("MB_SRG", "MB_SRS"),
              c(length(MB_SRG), length(MB_SRS)))
)
#Dice_group_set <- rbind(Dice_group_set, DiceB_group_set)
#Dice_group_set$group[Dice_group_set$group == "MW_SRG"] <- "Group"; Dice_group_set$group[Dice_group_set$group == "MB_SRG"] <- "Group"; Dice_group_set$group[Dice_group_set$group == "MW_SRS"] <- "Setting"; Dice_group_set$group[Dice_group_set$group == "MB_SRS"] <- "Setting"

Effects_group_set <- data.frame(
  value = c(EF_SRG, EF_SRS),
  group = rep(c("EF_SRG", "EF_SRS"),
              c(length(EF_SRG), length(EF_SRS)))
)

DGS <- ggplot(Dice_group_set, aes(x = value, fill = group)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 20) +
  theme_angele_ss+
  scale_x_continuous("Dice coefficients", breaks = seq(0.50, 0.60, by=0.01), limits = c(0.50, 0.60)) + 
  scale_y_continuous("N", breaks = seq(0, 50, by=10), limits = c(0, 50)) + 
  scale_fill_manual(labels = c("Groups", "Settings"), values = c("#007ABB", "#ba0f09")) +
  labs(fill = "Dice within:") +
  theme(legend.position = "bottom", legend.text = element_text(size = 12, vjust = 8), legend.title = element_text(size = 12, vjust = 0.6)) +
  ggtitle("(a)") 

DBGS <- ggplot(DiceB_group_set, aes(x = value, fill = group)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 20) +
  theme_angele_ss+
  scale_x_continuous("Dice coefficients", breaks = seq(0.45, 0.50, by=0.01), limits = c(0.45, 0.50)) + 
  scale_y_continuous("N", breaks = seq(0, 50, by=10), limits = c(0, 50)) + 
  scale_fill_manual(labels = c("Groups", "Settings"), values = c("#007ABB", "#ba0f09")) +
  labs(fill = "Dice within:") +
  theme(legend.position = "bottom", legend.text = element_text(size = 12, vjust = 8), legend.title = element_text(size = 12, vjust = 0.6)) +
  ggtitle("(a)") 

EFGS <- ggplot(Effects_group_set, aes(x = value, fill = group)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 20) +
  theme_angele_ss+
  scale_x_continuous("Effect values", breaks = seq(0.04, 0.12, by=0.01), limits = c(0.04, 0.12)) + 
  scale_y_continuous(" ", breaks = seq(0, 50, by=10), limits = c(0, 50)) + 
  scale_fill_manual(labels = c("Group", "Setting"), values = c("#007ABB", "#ba0f09")) +
  labs(fill = "Effect sizes for:") +
  theme(legend.position = "bottom", legend.text = element_text(size = 12, vjust = 8), legend.title = element_text(size = 12, vjust = 0.6)) +
  ggtitle("(b)") 

DiceSR <- ggarrange(DGS, EFGS)
DiceSR
