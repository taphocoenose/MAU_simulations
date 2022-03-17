### Skeletal parts and utility indices ###

# Code authored by Ryan Breslawski (rbreslawski@smu.edu)
# Last edited in R v 4.0.4 on a Windows 10 machine, Mar 17, 2022

# CODED FOR PARALLEL PROCESSING (detectCores()-1)

# load libraries
library(ggplot2)
library(patchwork)
library(parallel)
library(rethinking)

# Set seed for reproducibility
set.seed(98745)

# Create functions for simulation
source("SEAsim_SimFunctions.R")
# Create functions for plotting
source("SEAsim_PlotFunctions.R")

# Sequence of MNE values over which to simulate correlations
MNE_vals <- seq(5, 500, by=5)

# Number of simulated correlations per MNE value used to
# calculate exaggeration ratios and sign error rates. In
# the accompanying paper, this value is set to 100,000,
# allowing for high precision in estimates of sign error
# rates and exaggeration ratio distributions. However, 
# nCo values of this magnitude are recommended only if 
# users have access to a computing cluster that allows 
# programs to run for days or weeks. The default nCo 
# value of 10 allows for only rough estimates of sign 
# error rates and exaggeration ratio distributions, but 
# it is feasible to run the script on most personal 
# computers.
nCo <- 10

# Number of simulated correlations to plot per MNE value
n_perMNE <- 400

# Maximum MNE value to do nCo simulations. Any greater
# MNE value will only be simulated n_perMNE times for
# plotting purposes
MNE_sig_m <- 500

# vector of relationship names
rnames <- c("Strong","Moderate","Weak","None")

# Skeletal elements
elements <- c("mandible","atlas","axis","cervical vertebrae",
              "thoracic vertebrae","lumbar vertebrae","pelvis","ribs",
              "sternum","scapula","prox. humerus", "dist. humerus", 
              "prox. radio-cubitus", "dist. radio-cubitus", "carpals",
              "prox. metacarpal", "dist. metacarpal", "prox. femur",
              "dist. femur", "prox. tibia", "dist. tibia", "astragalus",
              "calcaneus", "prox. metatarsal", "dist. metatarsal", "phalanges")

# Anatomical frequencies (Lyman 1994:230). Carpals have been reduced from 12 to
# 2 to reflect sets of carpals rather than individual bones.
# 
#  Lyman, R.L. (1994). Vertebrate Taphonomy. Cambridge: Cambridge University 
#  Press.
freq <- c(2,1,1,5,13,6,1,26,7,rep(2,16),24)

# Standardized modified general utility index for caribou (Binford 1978:74). 
# 
#  Binford, L.R. (1978). Nunamiut Ethnoarchaeology. 2012 reprint. Clinton 
#  Corners: Eliot Werner. [Table 2.7]
utility <- c(13.89, 9.79, 9.79, 35.71, 45.53, 32.05, 47.89, 49.77, 64.13,
             43.47, 43.47, 36.52, 26.62, 22.23, 15.53, 12.18, 10.50, 100.00,
             100.00, 64.73, 47.09, 31.66, 31.66, 29.93, 23.93, 13.72)

# Create data frame to store element information for all 26 parts; 
# To be exported at the end of the script.
ElementsALL <- data.frame(Element=elements, Anatomical.Frequency=freq, 
                          Utility=utility, stringsAsFactors=FALSE)
# Append probabilities of recovery for each element, for each relationship
ElementsALL$Strong <- c(0.0093, 0.0012, 3e-04, 0.0556, 0.1526, 0.0403, 0.0167, 
                        0.2034, 0.1422, 0.0338, 0.0144, 0.0159, 0.0147, 0.0133, 
                        0.0029, 0.0075, 0.0019, 0.0688, 0.0399, 0.0222, 0.0223, 
                        0.0071, 0.0164, 0.0211, 0.0086, 0.0676)
ElementsALL$Moderate <- c(0.0073, 0.0021, 0.0032, 0.0219, 0.1533, 0.0596, 0.0146, 
                          0.3066, 0.0247, 0.0318, 0.0138, 0.0140, 0.0139, 0.0153, 
                          0.0112, 0.0080, 0.0080, 0.0604, 0.0103, 0.0486, 0.0078, 
                          0.0101, 0.0056, 0.0203, 0.0157, 0.1119)
ElementsALL$Weak <- c(0.0131, 0.0045, 0.0027, 0.0743, 0.2493, 0.0830, 0.0094, 
                      0.0920, 0.0798, 0.0173, 0.0218, 0.0286, 0.0258, 0.0170, 
                      0.0131, 0.0010, 0.0019, 0.0229, 0.0164, 0.009, 0.0184, 
                      0.0326, 0.0292, 0.0300, 0.0189, 0.0880)
ElementsALL$None <- freq/sum(freq)

# Offsets for Poisson regression
ElementsALL$os <- log(ElementsALL$Anatomical.Frequency)


# Long bones
LBelements <- c("humerus", "radio-cubitus", "metacarpal", "femur", 
                "tibia", "metatarsal")
# Standardized modified general utility index for caribou  long
# bones(Binford 1978:74). Values are rescaled to 0-100, since adding
# proximal and distal ends makes the scale 0-200.
# 
#  Binford, L.R. (1978). Nunamiut Ethnoarchaeology. 2012 reprint. Clinton 
#  Corners: Eliot Werner. [Table 2.7]
LBUtility <- c(43.47+36.52, 26.64+22.23, 12.18+10.50, 
               100.00+100.00, 64.73+47.09, 29.93+23.93)/2

# Create data frame to store element information for all 6 long bones; 
# To be exported at the end of the script.
ElementsLB <- data.frame(Element=LBelements, Anatomical.Frequency=rep(2, 6), 
                         Utility=LBUtility, stringsAsFactors=FALSE)


# Append probabilities of recovery for each element, for each relationship
ElementsLB$Strong <- c(0.1350, 0.0166, 0.0266, 0.3465, 0.3402, 0.1351)
ElementsLB$Moderate <- c(0.1721, 0.1797, 0.0965, 0.3007, 0.1763, 0.0747)
ElementsLB$Weak <- c(0.1787, 0.2171, 0.0599, 0.1870, 0.2013, 0.1560)
ElementsLB$None <- rep(2, 6)/6

# Offsets for Poisson regression
ElementsLB$os <- log(ElementsLB$Anatomical.Frequency)

PITable2 <- ElementsALL
PITable3 <- ElementsLB

for(i in 4:7){
  PITable2[,i] <- (PITable2[,i]/PITable2$Anatomical.Frequency)
  PITable2[,i] <- 100*(PITable2[,i]/max(PITable2[,i]))
  PITable3[,i] <- (PITable3[,i]/PITable3$Anatomical.Frequency)
  PITable3[,i] <- 100*(PITable3[,i]/max(PITable3[,i]))
}

# Save data for tables
write.csv(PITable2, "Table 2.csv", row.names=FALSE)
write.csv(PITable3, "Table 3.csv", row.names=FALSE)

# Hypothetical analyst table. Utility is Binford's (1978:74) 
# Modified General Utility Index for caribou.
PITable1 <- data.frame(part=c("Humerus", "Radius", "Metacarpal",
                              "Femur", "Tibia", "Metatarsal", "TMNE", 
                              "rho critical value", "Sample rho"),
                       Utility=c(43.47+36.52, 26.64+22.23, 12.18+10.50,
                                 100+100, 64.73+47.09, 29.93+23.93, NA,
                                 NA, NA),
                       MNE1=c(1, 2, 0, 4, 2, 1, 10, 0.886, NA),
                       MNE2=c(3, 2, 4, 6, 4, 1, 20, 0.886, NA),
                       MNE4=c(7, 4, 5, 14, 8, 2, 40, 0.886, NA),
                       stringsAsFactors=FALSE)
PITable1[9,3:5] <- lapply(3:5, function(x){
  sresult <- cor.test(PITable1[1:6,x]/2, PITable1[1:6,2], 
           method="spearman", exact=FALSE)
  return(sresult$estimate)
})
write.csv(PITable1, "Table 1.csv", row.names=FALSE)

# Obtain population plots and statistics
popALL <- PopBones(ElementsALL, "26 parts", c("a", "c", "e"))
popLB <- PopBones(ElementsLB, "6 long bones", c("b", "d", "f"))

# Export plots of population relationships
PopPlot <- popALL$PopPlots[[1]] + popLB$PopPlots[[1]] +
  popALL$PopPlots[[2]] + popLB$PopPlots[[2]] +
  popALL$PopPlots[[3]] + popLB$PopPlots[[3]] +
  plot_layout(ncol=2)
# Suppress warning: ends of regression cut off in panel C
suppressWarnings(ggsave("Figure 1.jpeg", plot=PopPlot , 
                        device="jpeg", units="in", width=9, 
                        height=7, dpi=600))

# Perform simulations for each set of elements
simALL <- SimBones(x1=ElementsALL, cDF=popALL$cDF, MNE_values=MNE_vals,
                   nCor=nCo, n_per_MNE=n_perMNE, 
                   MNE_sig_max=MNE_sig_m, rnames1=rnames)
simLB <- SimBones(x1=ElementsLB, cDF=popLB$cDF, MNE_values=MNE_vals,
                  nCor=nCo, n_per_MNE=n_perMNE, 
                  MNE_sig_max=MNE_sig_m, rnames1=rnames)

# Save output to working directory
save(file="Results.RData", simALL, simLB, popALL, 
     popLB, ElementsALL, ElementsLB, rnames)


# Create plots for simulated values
pALLrho <- SimPlot(simALL$SimsValues, "26 parts", popALL$cDF,
                   "rho", c("a", "c", "e", "g"),  rnames, T)
pLBrho <- SimPlot(simLB$SimsValues, "6 long bones", popLB$cDF,
                  "rho", c("b", "d", "f", "h"), rnames, F)
pALLbeta <- SimPlot(simALL$SimsValues, "26 parts", popALL$cDF,
                    "beta", c("a", "c", "e", "g"), rnames, T)
pLBbeta <- SimPlot(simLB$SimsValues, "6 long bones", popLB$cDF,
                   "beta", c("b", "d", "f", "h"), rnames, F)
pALLBayes <- SimPlot(simALL$SimsValues, "26 parts", popALL$cDF,
                     "beta_Bayes", c("a", "c", "e", "g"), rnames, T)
pLBBayes <- SimPlot(simLB$SimsValues, "6 long bones", popLB$cDF,
                    "beta_Bayes", c("b", "d", "f", "h"), rnames, F)

# Panel layout for simulated correlation plot panels
layout1 <- "
AB
CD
EF
GH
II"

plot_rho <- pALLrho[[1]] + pLBrho[[1]] + pALLrho[[2]] + 
  pLBrho[[2]] + pALLrho[[3]] + pLBrho[[3]] + pALLrho[[4]] + 
  pLBrho[[4]] + pALLrho[[5]] + 
  plot_layout(design=layout1, heights=c(rep(4, 4), 1))
ggsave("Figure 2.jpeg", plot=plot_rho, device="jpeg", units="in",
       width=9, height=9, dpi=600)

plot_beta <- pALLbeta[[1]] + pLBbeta[[1]] + pALLbeta[[2]] + 
  pLBbeta[[2]] + pALLbeta[[3]] + pLBbeta[[3]] + pALLbeta[[4]] + 
  pLBbeta[[4]] + pALLbeta[[5]] +
  plot_layout(design=layout1, heights=c(rep(4, 4), 1))
ggsave("Figure 5.jpeg", plot=plot_beta, device="jpeg", units="in",
       width=9, height=9, dpi=600)

plot_Bayes <- pALLBayes[[1]] + pLBBayes[[1]] + pALLBayes[[2]] + 
  pLBBayes[[2]] + pALLBayes[[3]] + pLBBayes[[3]] + pALLBayes[[4]] + 
  pLBBayes[[4]] +  
  plot_layout(design=layout1, heights=c(rep(4, 4), 1))
ggsave("Figure 8.jpeg", plot=plot_Bayes, device="jpeg", units="in",
       width=9, height=9, dpi=600)

# Create exaggeration plots
pexALLrho <- ExagPlot(simALL$SimsSummary, "26 parts", popALL$cDF,  
                      "rho", letters[1:3], rnames[1:3], F, c(0, 3))
pexLBrho <- ExagPlot(simLB$SimsSummary, "6 long bones", popLB$cDF,  
                     "rho", letters[4:6], rnames[1:3], T, c(0, 3))
pexALLbeta <- ExagPlot(simALL$SimsSummary, "26 parts", popALL$cDF,  
                       "beta", letters[1:3], rnames[1:3], F, c(0, 7))
pexLBbeta <- ExagPlot(simLB$SimsSummary, "6 long bones", popLB$cDF,  
                      "beta", letters[4:6], rnames[1:3], T, c(0, 7))

layout2 <- "
ABC
DEF
MMM"

plot_exag <-  pexALLrho[[1]] + pexALLrho[[2]] + pexALLrho[[3]] + 
  pexLBrho[[1]] + pexLBrho[[2]] + pexLBrho[[3]] + pexALLrho[[4]] +
  plot_layout(design=layout2, heights=c(4, 4, 1))
ggsave("Figure 3.jpeg", plot=plot_exag, device="jpeg", units="in",
       width=8, height=5, dpi=600)


plot_exag_beta <- pexALLbeta[[1]] + pexALLbeta[[2]] + 
  pexALLbeta[[3]] + pexLBbeta[[1]] + pexLBbeta[[2]] +
  pexLBbeta[[3]] + pexALLbeta[[4]] +
  plot_layout(design=layout2, heights=c(rep(4, 2), 1))
ggsave("Figure 6.jpeg", plot=plot_exag_beta, device="jpeg",
       units="in", width=8, height=5, dpi=600)


# Create sign error plots
psignALL <- SignPlot(simALL$SimsSummary, "26 parts", popALL$cDF, 
                     rep("a", 3), c(T, T, T), c(T, T, T), 0.4)
psignLB <- SignPlot(simLB$SimsSummary, "6 long bones", popLB$cDF, 
                    rep("b", 3), c(T, T, T), c(F, F, F), 0.4)

plot_sign <- psignALL[[2]] +  psignLB[[2]] + plot_layout(ncol=2)
ggsave("Figure 4.jpeg", plot=plot_sign, device="jpeg",
       units="in", width=9, height=3, dpi=600)

plot_sign_beta <- psignALL[[3]] + psignLB[[3]] +  plot_layout(ncol=2)
ggsave("Figure 7.jpeg", plot=plot_sign_beta, device="jpeg",
       units="in", width=9, height=3, dpi=600)


# Create plots for samples of simulated Bayesian results
pbayesALL <- BayesPlot(simALL$SampleBayes, simLB$SampleBayes,
                       popALL$cDF, "26 parts", 500, 50,
                       rnames, c(T, T, T, T), c(F, F, F, T),
                       c("a", "c", "e", "g"))

pbayesLB <- BayesPlot(simLB$SampleBayes, simALL$SampleBayes,
                      popLB$cDF, "6 long bones", 500, 50, 
                      rnames, c(F, F, F, F), c(F, F, F, T),
                      c("b", "d", "f", "h"))

plot_BSamples <- pbayesALL[[1]] + pbayesLB[[1]] + pbayesALL[[2]] + 
  pbayesLB[[2]] + pbayesALL[[3]] + pbayesLB[[3]] + pbayesALL[[4]] + 
  pbayesLB[[4]] +  
  plot_layout(design=layout1, heights=c(rep(4, 4), 1))
ggsave("Figure 9.jpeg", plot=plot_BSamples, device="jpeg", units="in",
       width=9, height=9, dpi=600)

# Subset simulation summary data for ease with text description
# of results
subMNE <- c(5, 20, 50, 100, 250, 500)
TextDescALL <- simALL$SimsSummary[which(simALL$SimsSummary$MNE%in%subMNE),]
TextDescLB <- simLB$SimsSummary[which(simLB$SimsSummary$MNE%in%subMNE),]
TextDesc <- rbind(TextDescALL, TextDescLB)
rm(TextDescALL, TextDescLB)
TextDesc <- TextDesc[-which(TextDesc$Relationship=="None"),]
TextDesc$sign_rho <- round(TextDesc$sign_rho, 3)
TextDesc$sign_sig_rho <- round(TextDesc$sign_sig_rho, 3)
TextDesc$sign_beta <- round(TextDesc$sign_beta, 3)
TextDesc$sign_sig_beta <- round(TextDesc$sign_sig_beta, 3)
TextDesc$Elements <- c(rep("All", nrow(TextDesc)/2), 
                        rep("LB", nrow(TextDesc)/2))
TextDesc$rho <- paste0(round(TextDesc$rho_med, 3), " (", 
                       round(TextDesc$rho_lo, 3), "-",
                       round(TextDesc$rho_up, 3), ")")
TextDesc$rho_sig <- paste0(round(TextDesc$sig_rho_med, 3), " (", 
                       round(TextDesc$sig_rho_lo, 3), "-",
                       round(TextDesc$sig_rho_up, 3), ")")
TextDesc$beta <- paste0(round(TextDesc$beta_med, 3), " (", 
                       round(TextDesc$beta_lo, 3), "-",
                       round(TextDesc$beta_up, 3), ")")
TextDesc$beta_sig <- paste0(round(TextDesc$sig_beta_med, 3), " (", 
                           round(TextDesc$sig_beta_lo, 3), "-",
                           round(TextDesc$sig_beta_up, 3), ")")
ResultsTable <- TextDesc[,c(1, 2, 33:37, 28, 31, 29, 32)]

write.csv(ResultsTable, "Tables 4-7.csv", row.names=F)
