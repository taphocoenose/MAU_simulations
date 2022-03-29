### Skeletal parts and utility indices ###

# Code authored by Ryan Breslawski (rbreslawski@smu.edu)
# Last edited in R v 4.0.4 on a Windows 10 machine, Mar 22, 2022

# load libraries
library(ggplot2)
library(patchwork)
library(rethinking)

# Set seed for reproducibility
set.seed(98745)

# Create functions for script
source("SEAexamples_Functions.R")

# Create data frames for the two archaeological examples. Manzanilla
# Site data are from Delsol and Gerouard 2016, and Tel Kabri data are
# from Marom et al. 2015. Mandibles are excluded for Tel Kabri,
# as there are no density values for these elements.
# 
#  Delsol, N., & Grouard, S. (2016). Comments on Amerindian hunting 
#  practices in Trinidad (West Indies): Tetrapods from the Manzanilla 
#  Site (Late Ceramic age 300-900 AD). The Journal of Island and Coastal 
#  Archaeology, 11, 385-410.
#  
#  Marom, N., Yasur-Landau, A., & Cline, E.H. (2015). The silent coast: 
#  Zooarchaeological evidence to the development of a second millennium 
#  palace at Tel Kabri. Journal of Anthropological Archaeology, 39, 
#  181-192.

E_Manz <- data.frame(Element=c("Femur", "Tibia", "Metatarsal", "Humerus", 
                               "Radius", "Mandible", "Skull", "Metacarpal"),
                     AF=c(2, 2, 2, 2, 2, 2, 1, 2),
                     cerv_MNE=c(19, 7, 3, 6, 10, 6, 4, 4),
                     pecc_MNE=c(5, 6, 1, 4, 3, 5, 7, 3),
                     SFUI=c(100, 62.8, 37, 36.8, 25.8, 11.5, 9.1, 5.2),
                     stringsAsFactors=FALSE)
E_TelK <- data.frame(Element=c("Cervical vert",
                               "Thoracic vert", "Lumbar vert", "Pelvis",
                               "Femur", "Tibia", "Metatarsus", "Scapula",
                               "Humerus", "Radius", "Metacarpus", 
                               "Astragalus", "Calcaneus", "Phalanx 1",
                               "Phalanx 2", "Phalanx 3"),
                     AF=c(7, 13, 7, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                8, 8, 8),
                     EP_MNE=c(4, 2, 4, 5, 5, 7, 4, 2, 15, 6, 4, 3, 3, 
                              12, 3, 3),
                     LP_MNE=c(14, 15, 9, 7, 13, 10, 12, 9, 15, 16, 13, 
                              10, 8, 48, 17, 10),
                     Meat=c(55.32, 46.47, 38.88, 81.3, 78.24, 20.76,
                            6.37, 44.89, 28.24, 14.01, 4.74, 6.37, 6.37, 
                            3.37, 3.37, 3.37),
                     Marrow=c(1, 1, 1, 9.57, 47, 79, 69, 6.23, 34, 52,
                              67.5, 1, 23.11, 33.77, 25.11, 1),
                     Density=c(0.11, 0.24, 0.22, 0.26, 0.36, 0.59, 0.53,
                               0.25, 0.37, 0.52, 0.55, 0.54, 0.56, 0.4, 
                               0.39, 0.3),
                     stringsAsFactors=FALSE)


# Put density values on 1-100 scale
E_TelK$Density <- 100*(E_TelK$Density/max(E_TelK$Density))

# Create data frame for data from Morin et al. (2017a,b).
# 
#  Morin, E.A., Ready, E., Boileau, A., Beaval, C., & Coumont, M.-P. 
#  (2017a). Problems of identification and quantification in a
#  archaeozoological analysis, part I: insights from a blind test.
#  Journal of Archaeological Method and Theory, 24, 886-937.
#  
#  #  Morin, E.A., Ready, E., Boileau, A., Beaval, C., & Coumont, M.-P. 
#  (2017a). Problems of identification and quantification in a
#  archaeozoological analysis, part II: presentation of an alternative
#  souting method. Journal of Archaeological Method and Theory, 24, 
#  886-937.

MorinEtAl <- data.frame(Element=c("humerus", "radius", "ulna",
                                  "metacarpal", "femur", "tibia",
                                  "metatarsal", "cranium", "mandible",
                                  "rib", "scapula", "innominate",
                                  "malleolus", "talus", "calcaneus",
                                  "phalanx 1", "phalanx 2", "phalanx 3"),
                        AF=c(2, 2, 2, 2, 2, 2, 2, 1, 2, 26, 2, 2, 2, 2,
                             2, 8, 8, 8),
                        MCE_ANE=c(23, 21, 21, 26, 24, 21, 26, 1, 2, 2,
                                  7, 3, 7, 3, 3, 2, 2, 3),
                        BGRE_ANE=c(20, 20, 20, 24, 29, 21, 21, 1, 3, 2,
                                   6, 2, 7, 2, 3, 1, 3, 2),
                        A_MCE=c(23, 22, 22, 26, 22, 20, 25, 1, 2, 2, 5, 
                                2, 6, 2, 3, 2, NA, 2),
                        B_MCE=c(23, 21, 21, 25, 25, 20, 22, 3, 2, 2, 4, 
                                3, 7, 2, 3, 2, 2, 2),
                        D_MCE=c(23, 21, 22, 23, 22, 23, 23, 1, 2, 2, 4,
                                2, 7, 2, 2, 1, 1, 0),
                        NDE_MCE=c(23, 21, 19, 26, 20, 19, 22, NA, 2, 2,
                                  4, 2, 6, 2, 2, 2, 1, 2),
                        A_BGRE=c(21, 20, 20, 22, 29, 22, 22, 1, 3, 1, 
                                 3, 2, 6, 2, 3, 1, 2, 2),
                        B_BGRE=c(27, 20, 22, 23, 31, 23, 24, NA, 2, 1,
                                 2, 2, NA, NA, NA, 2, 2, 2),
                        D_BGRE=c(22, 17, 21, 24, 29, 22, 19, 1, 3, 1,
                                 3, 2, 5, 3, 3, 2, 2, 4),
                        NDE_BGRE=c(17, 15, 16, 19, 20, 18, 18, NA, 2, 1,
                                   1, 2, 4, 2, 2, 1, 3, 1),
                        stringsAsFactors=FALSE)

# Create index predictor
MorinEtAlUtilityIndex <- c(65, 70, 85, 90, 60, 51, 95, 14, 10, 
                           18, 45, 20, 31, 15, 25, 3, 7, 1)
MorinEtAl$Analyst_A_MNE <- MorinEtAl$Known_sample <- MorinEtAlUtilityIndex
MorinEtAl$Analyst_D_MNE <- MorinEtAl$Analyst_B_MNE <- MorinEtAl$Analyst_A_MNE
MorinEtAl$NDE <- MorinEtAl$Analyst_D_MNE
  
# Subset to remove missing value cases
MEA_MC <- MorinEtAl[!(is.na(MorinEtAl$A_MCE) | 
                        is.na(MorinEtAl$NDE_MCE)),]
MEA_BGR <- MorinEtAl[!(is.na(MorinEtAl$B_BGRE) | 
                         is.na(MorinEtAl$NDE_BGRE)),]

write.csv(MorinEtAl, "Table 9.csv", row.names=FALSE)

# Fit models
E_Manz_pecc <- BayesModel(E_Manz, "pecc_MNE", "SFUI", "AF")
E_Manz_cerv <- BayesModel(E_Manz, "cerv_MNE", "SFUI", "AF")
E_TelK_EP_dens <- BayesModel(E_TelK, "EP_MNE", "Density", "AF")
E_TelK_EP_meat <- BayesModel(E_TelK, "EP_MNE", "Meat", "AF")
E_TelK_EP_marrow <- BayesModel(E_TelK, "EP_MNE", "Marrow", "AF")
E_TelK_EP_densmeat <- BayesModel(E_TelK, "EP_MNE", c("Density", "Meat"), "AF")
E_TelK_EP_densmar <- BayesModel(E_TelK, "EP_MNE", c("Density", "Marrow"), "AF")
E_TelK_LP_dens <- BayesModel(E_TelK, "LP_MNE", "Density", "AF")
E_TelK_LP_meat <- BayesModel(E_TelK, "LP_MNE", "Meat", "AF")
E_TelK_LP_marrow <- BayesModel(E_TelK, "LP_MNE", "Marrow", "AF")
E_TelK_LP_densmeat <- BayesModel(E_TelK, "LP_MNE", c("Density", "Meat"), "AF")
E_TelK_LP_densmar <- BayesModel(E_TelK, "LP_MNE", c("Density", "Marrow"), "AF")
MorinEtAl_MC_ANE <- BayesModel(MEA_MC, "MCE_ANE", "Known_sample", "AF")
MorinEtAl_MC_A <- BayesModel(MEA_MC, "A_MCE", "Analyst_A_MNE", "AF")
MorinEtAl_MC_B <- BayesModel(MEA_MC, "B_MCE", "Analyst_B_MNE", "AF")
MorinEtAl_MC_D <- BayesModel(MEA_MC, "D_MCE", "Analyst_D_MNE", "AF")
MorinEtAl_MC_NDE <- BayesModel(MEA_MC, "NDE_MCE", "NDE", "AF")
MorinEtAl_BGR_ANE <- BayesModel(MEA_BGR, "BGRE_ANE", "Known_sample", "AF")
MorinEtAl_BGR_A <- BayesModel(MEA_BGR, "A_BGRE", "Analyst_A_MNE", "AF")
MorinEtAl_BGR_B <- BayesModel(MEA_BGR, "B_BGRE", "Analyst_B_MNE", "AF")
MorinEtAl_BGR_D <- BayesModel(MEA_BGR, "D_BGRE", "Analyst_D_MNE", "AF")
MorinEtAl_BGR_NDE <- BayesModel(MEA_BGR, "NDE_BGRE", "NDE", "AF")

# Extract samples from Manzanilla posteriors to calculate contrast
Man_pecc_B <- extract.samples(E_Manz_pecc$model)$bu
Man_cerv_B <- extract.samples(E_Manz_cerv$model)$bu
Man_con <- quantile(Man_cerv_B - Man_pecc_B, 
                    probs=c(0.5, 0.025, 0.975))

# Create plots for posterior parameters
plot_pars_M <- PostParsPlot(list(E_Manz_pecc, E_Manz_cerv), 
                            c("springgreen4", "purple1"), 
                            "", TRUE, "Model", "right",
                            c("Peccary", "Cervid"))
plot_pars_TEP <- PostParsPlot(list(E_TelK_EP_dens, E_TelK_EP_meat, 
                                   E_TelK_EP_marrow, E_TelK_EP_densmeat, 
                                   E_TelK_EP_densmar), 
                               c("grey64", "red", "blue",
                                 "tomato2", "deepskyblue2"), 
                              "a", FALSE, "Model", "right")
plot_pars_TLP <- PostParsPlot(list(E_TelK_LP_dens, E_TelK_LP_meat, 
                                   E_TelK_LP_marrow, E_TelK_LP_densmeat, 
                                   E_TelK_LP_densmar), 
                              c("grey64", "red", "blue",
                                "tomato2", "deepskyblue2"), 
                              "b", TRUE, "", "none")
plot_pars_MCE <- PostParsPlot(list(MorinEtAl_MC_ANE, MorinEtAl_MC_D, 
                                   MorinEtAl_MC_B, MorinEtAl_MC_A,
                                   MorinEtAl_MC_NDE), 
                              rep("grey64", 5), "a", TRUE, "Analyst", "none")
plot_pars_BGRE <- PostParsPlot(list(MorinEtAl_BGR_ANE, MorinEtAl_BGR_D, 
                                   MorinEtAl_BGR_B, MorinEtAl_BGR_A,
                                   MorinEtAl_BGR_NDE), 
                              rep("grey64", 5), "b", TRUE, "Analyst", "none")

ggsave("Figure 11.jpeg", plot_pars_M, "jpeg", width=7, height=2,
       units="in", dpi=600)

TK_pars_plot <- plot_pars_TEP + plot_pars_TLP + plot_layout(heights=c(2, 2))

ggsave("Figure 13.jpeg", TK_pars_plot, "jpeg", width=7, height=4,
       units="in", dpi=600)


Analyst_comparison_plot <- plot_pars_MCE + plot_pars_BGRE + plot_layout(nrow=1)

ggsave("Figure 15.jpeg", Analyst_comparison_plot, "jpeg", width=11, height=4,
       units="in", dpi=600)

# Create plots for relationships
Rel_plots <- list(Manz_pecc=PostModelPlot(E_Manz_pecc, "a"), 
                  Manz_cerv=PostModelPlot(E_Manz_cerv, "b"),
                  TelKEP_dens=PostModelPlot(E_TelK_EP_dens, "a"),
                  TelKLP_dens=PostModelPlot(E_TelK_LP_dens, "b"),
                  TelKEP_meat=PostModelPlot(E_TelK_EP_meat, "c"),
                  TelKLP_meat=PostModelPlot(E_TelK_LP_meat, "d"),
                  TelKEP_marrow=PostModelPlot(E_TelK_EP_marrow, "e"),
                  TelKLP_marrow=PostModelPlot(E_TelK_LP_marrow, "f"))

Rel_plot_Man <- Rel_plots[[1]] + Rel_plots[[2]] 
Rel_plot_Tel <- Rel_plots[[3]] + Rel_plots[[4]] + Rel_plots[[5]] + 
                Rel_plots[[6]] + Rel_plots[[7]] + Rel_plots[[8]] + 
                plot_layout(ncol=2)

ggsave("Figure 12.jpeg", Rel_plot_Man, "jpeg", width=9, height=3,
       units="in", dpi=600)
ggsave("Figure 14.jpeg", Rel_plot_Tel, "jpeg", width=9, height=9,
       units="in", dpi=600)

# Model comparisons
compare_TEP <- compare(E_TelK_EP_dens$model, E_TelK_EP_densmeat$model, 
                       E_TelK_EP_densmar$model, E_TelK_EP_meat$model, 
                       E_TelK_EP_marrow$model)
compare_TLP <- compare(E_TelK_LP_dens$model, E_TelK_LP_densmeat$model, 
                       E_TelK_LP_densmar$model, E_TelK_LP_meat$model, 
                       E_TelK_LP_marrow$model)
suppressWarnings(write.csv(compare_TEP, "Table 8 [part 2].csv"))
suppressWarnings(write.csv(compare_TLP, "Table 8 [part 3].csv"))

# Extract parameter summaries from all models and combine them
allmodels <- list(E_Manz_cerv, E_Manz_pecc, E_TelK_EP_dens, 
                  E_TelK_EP_densmar, E_TelK_EP_densmeat, 
                  E_TelK_EP_marrow, E_TelK_EP_meat, 
                  E_TelK_LP_dens, E_TelK_LP_densmar, 
                  E_TelK_LP_densmeat, E_TelK_LP_marrow, 
                  E_TelK_LP_meat)

allmodelsDF <- lapply(allmodels, function(x){
  mname <- paste(x$MNE_vals, paste(x$predictor, collapse=" "))
  df <- x$summary
  df$model <- rep(mname, nrow(df))
  return(df)
})

allmodelsDF <- do.call("rbind", allmodelsDF)

suppressWarnings(write.csv(allmodelsDF, "Table 8 [part 1].csv", 
                           row.names=FALSE))

# Create plots illustrating hypothetical priors, likelihoods,
# and posteriors in two data scenarios

# x-grid values to plot distributions
grid <- seq(-0.07, 0.1, length.out=200)

prior <- dnorm(grid, 0, 0.01)
prior <- prior/sum(prior)

likelihood <- matrix(ncol=2,
                     c(dnorm(grid, 0.04, 0.016),
                     dnorm(grid, 0.024, 0.005)))
posterior <- sapply(1:2, function(x){
  ustpost <- likelihood[,x]*prior
  return(ustpost/sum(ustpost))
})

prior <- prior/max(posterior)
posterior <- posterior/max(posterior)
likelihood <- max(posterior)*(likelihood/max(likelihood))

p <- lapply(1:2, function(i){
  
  plt <- ggplot(NULL)+
    annotate("text", label=letters[i], x=-0.06,
             y=max(posterior)*0.9)+
    annotate("segment", x=0.02, xend=0.02, y=Inf, yend=-Inf,
             color="grey")+
    annotate("line", x=grid, y=prior, color="blue",
             linetype="dotted")+
    annotate("line", x=grid, y=likelihood[,i], color="red",
             linetype="dashed")+
    annotate("line", x=grid, y=posterior[,i], color="purple")+
    scale_y_continuous(limits=c(0,1))+
    labs(x="Parameter value")+
    theme(panel.grid=element_blank(), 
          axis.text.y=element_blank(),
          panel.background=element_blank(),
          panel.border=element_rect(fill=NA, color="grey"),
          axis.ticks.x=element_line(color="grey"),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())
  
  return(plt)
})

ggsave("Figure 8.jpeg", p[[1]] + p[[2]], "jpeg", 
       width=7, height=2, units="in", dpi=600)


save.image(file="Results_examples.RData")

