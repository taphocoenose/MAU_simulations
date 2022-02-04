### Skeletal parts and utility indices ###

# Code authored by Ryan Breslawski (rbreslawski@smu.edu)
# Last edited in R v 4.0.4 on a Windows 10 machine, Jan 25, 2022

# These functions are required for and called within SEAsim.R. They rely
# on functions from the ggplot2 package, which is called in SEAsim.R

# Function to plot raw simulated values
SimPlot <- function(sDF, bones, truevalDF, cortype, plotlets, 
                    plotRels, ylabson){
  
  # Set colors for points
  if(cortype != "beta_Bayes"){
    sc <- c("orange", "grey")
  }else{
    sc <- rep("orchid3", 2)
  }
  
  # x max for plotting
  MNE_max <- max(sDF$MNE)*1.01
  
  # Find correlation for plotting
  sDF$corVal <- sDF[, which(colnames(sDF)==cortype)]
  
  if(cortype!= "beta_Bayes"){
    # Find p value for plotting
    sDF$p <- sDF[, which(colnames(sDF)==paste0(cortype, "_p"))]
    cortype1 <- cortype
    if(cortype != "beta"){
      yinc <- seq(1, -1, by=-0.5)
    }else{
      yinc <- seq(0.05, -0.05, by=-0.025)
    }
  }else{
    # Dummy p value column for plotting
    sDF$p <- sample(c("y", "n"), nrow(sDF), replace=TRUE)
    cortype1 <- "beta"
    yinc <- seq(0.036, -0.036, by=-0.018)
  }
  
  simplots <- lapply(plotRels, function(x){
    
    sDF <- sDF[which(sDF$Relationship==x), ]
    
    # If beta_Bayes, get median value for each MNE
    if(cortype=="beta_Bayes"){
      bmed <- sapply(unique(sDF$MNE), function(y){
        c(median(sDF$corVal[which(sDF$MNE==y)]), y)
      })
    }
    
    # If any outlying corVal values are beyond the y boundaries
    # of the plot, remove them
    sDF <- sDF[which(sDF$corVal <= max(yinc)),]
    sDF <- sDF[which(sDF$corVal >= min(yinc)),]
    
    # Subset true value from truevalDF and format for plotting
    if(x != "None"){
      truevalue <- truevalDF[which(truevalDF$str==x),
                             which(colnames(truevalDF)==cortype1)]
      truevalString <- as.character(round(truevalue, 3))
      if(nchar(truevalString < 5)){
        zeros <- paste0(rep("0", (5-nchar(truevalString))), collapse="")
        truevalString <- paste0(truevalString, zeros)
      }
    }else{
      truevalue <- 0
      truevalString <- "0.000"
    }
    
    if(bones=="26 parts"){
      bones2 <- "26 part types"
    }else{
      bones2 <- "6 part types"
    }
  
    if(cortype%in%c("beta", "beta_Bayes")){
      cortypeex <- parse(text="expression(italic(beta))")
      truevalString <- paste0("italic(beta)=='", 
                             truevalString, "'")
    }else if(cortype=="rho"){
      cortypeex <- parse(text="expression(italic('rho'))")
      truevalString <- paste0("italic('rho')=='", 
                             truevalString, "'")
    }else{
      cortypeex <- parse(text="expression(italic('r'))")
      truevalString <- paste0("italic('r')=='", 
                             truevalString, "'")
    }
    
    # Create plot base
    simplot <- ggplot(sDF, aes(x=MNE, y=corVal, color=p))+
      annotate("segment",x=-Inf, xend=Inf, y=0, yend=0, size=0.5)+
      geom_point(alpha=0.1, shape=16, size=1)
    
    # If plotting Bayesian relationships, include median
    # value for beta
    if(cortype=="beta_Bayes"){
      simplot <- simplot+
        annotate("line", x=bmed[2,], y=bmed[1,],
                 size=0.7, color="white")
    }
    
    xtitle <- ifelse(x==plotRels[length(plotRels)], 
                     parse(text="element_text()"),
                     parse(text="element_blank()"))
    ytitle <- ifelse(ylabson, 
                     parse(text="element_text()"),
                     parse(text="element_blank()"))
    
    # Add final plot elements
    simplot <- simplot+
      annotate("segment",x=-Inf, xend=Inf, y=truevalue, 
               yend=truevalue, size=0.7, linetype=2)+
      scale_color_manual(values=c("y"=sc[1], "n"=sc[2]))+
      scale_x_continuous(limits=c(0, MNE_max),
                       expand=c(0.0175,0.0175))+ 
      scale_y_continuous(limits=1.1*c(min(yinc), max(yinc)),
                         breaks=yinc, expand=c(0,0))+
      labs(x="TMNE", y=eval(cortypeex))+
      annotate("text",label=plotlets[which(plotRels==x)], 
               x=MNE_max, y=min(yinc), hjust=1, vjust=0, size=6)+
      annotate("text",label=truevalString, x=MNE_max*0.75, 
               y=min(yinc)+0.002*(max(yinc)-min(yinc)), 
               vjust=0, size=4, parse=TRUE)+
      annotate("text",label=bones2, x=MNE_max*0.75, 
               y=min(yinc)+0.2*(max(yinc)-min(yinc)),
               vjust=1, size=4)+
      theme(panel.background=element_blank(),
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            axis.title.x=eval(xtitle),
            axis.title.y=eval(ytitle),
            panel.border=element_rect(color="black", 
                                      fill=NA, size=0.6),
            legend.position="none")
    
    return(simplot)
    
  })
  
  # Create plot legend for plots that show significant vs not
  # significant correlations
  if(cortype != "beta_Bayes"){
    # Simulate random points for key symbols.
    kdf <- data.frame(y=rbeta(500,2,2), xsig=2.7-rbeta(500,1,2),
                      xnsig=5-rbeta(500,1,2))
    # Create plot
    plotkey <- ggplot(data=kdf, aes(y=y))+
      geom_point(aes(x=xsig), color=sc[2], size=1, alpha=0.15)+
      geom_point(aes(x=xnsig), color=sc[1], size=1, alpha=0.1)+
      annotate("text", label="p > 0.05", color=sc[2], size=5,
               hjust=0, x=2.74, y=0.5)+
      annotate("text", label="p < 0.05", color=sc[1], size=5,
               hjust=0, x=5.04, y=0.5)+
      scale_x_continuous(limits=c(0,7),expand=c(0,0))+
      scale_y_continuous(limits=c(-1.5,2.5),expand=c(0,0))+
      theme_void()
    
    simplots[[length(simplots) + 1]] <- plotkey
  }
  
  return(simplots)

}

# Function to plot exaggeration ratios
ExagPlot <- function(sDF, bones, truevalDF, cortype, plotlets, 
                    plotRels, xlabson, ylims){
  
  # Remove simulations for which exaggeration ratios
  # were not calculated due to high MNE values or 
  # "None" relationship simulations
  sDF <- sDF[which(!is.na(sDF$r_med)), ]
  
  # Colors for significant and not significant plot symbols
  sc <- c("orange", "grey")
  
  # Find and relabel columns for plotting based on
  # cor type
  sDF$nsM <- sDF[, which(colnames(sDF)==paste0(cortype, "_med"))]
  sDF$nsU <- sDF[, which(colnames(sDF)==paste0(cortype, "_up"))]
  sDF$nsL <- sDF[, which(colnames(sDF)==paste0(cortype, "_lo"))]
  sDF$sM <- sDF[, which(colnames(sDF)==paste0("sig_", cortype, "_med"))]
  sDF$sU <- sDF[, which(colnames(sDF)==paste0("sig_", cortype, "_up"))]
  sDF$sL <- sDF[, which(colnames(sDF)==paste0("sig_", cortype, "_lo"))]
  
  explots <- lapply(plotRels, function(x){
    
    sDF <- sDF[which(sDF$Relationship==x), ]
    
    # If any outlying corVal values are beyond the ymax boundary
    # of the plot, truncate them
    sDF$nsU[which(sDF$nsU > ylims[2])] <- ylims[2]
    sDF$sU[which(sDF$sU > ylims[2])] <- ylims[2]
    sDF$nsM[which(sDF$nsM > ylims[2])] <- ylims[2]
    sDF$sM[which(sDF$sM > ylims[2])] <- ylims[2]
    
    # Subset true value from truevalDF and format for plotting
    truevalue <- truevalDF[which(truevalDF$str==x),
                           which(colnames(truevalDF)==cortype)]
    truevalString <- as.character(round(truevalue, 3))
    if(nchar(truevalString < 5)){
      zeros <- paste0(rep("0", (5-nchar(truevalString))), collapse="")
      truevalString <- paste0(truevalString, zeros)
    }
    
    MNE_max <- max(sDF$MNE)
    
    xtitle <- ifelse(xlabson, 
                     parse(text="element_text()"),
                     parse(text="element_blank()"))
    ytitle <- ifelse(x==plotRels[1], 
                     parse(text="element_text()"),
                     parse(text="element_blank()"))
    
    if(bones=="26 parts"){
      bones2 <- "26 part types"
    }else{
      bones2 <- "6 part types"
    }
    
    if(cortype%in%c("beta", "beta_bayes")){
      truevalString <- paste0("italic(beta)=='", truevalString, "'")
    }else if(cortype=="rho"){
      truevalString <- paste0("italic('rho')=='",  truevalString, "'")
    }else{
      truevalString <- paste0("italic('r')=='", truevalString, "'")
    }
    
    xtext <- (MNE_max - min(sDF$MNE))*0.05 + min(sDF$MNE)
    
    # Create plot base
    explot <- ggplot(sDF, aes(x=MNE))+
      annotate("segment", x=-Inf, xend=Inf, y=1, yend=1, 
               size=0.8, linetype=2)+
      geom_ribbon(aes(ymin=nsL, ymax=nsU), fill=sc[2], alpha=0.2)+
      geom_line(aes(y=nsM), color=sc[2], size=0.8)+
      geom_ribbon(aes(ymin=sL, ymax=sU), fill=sc[1], alpha=0.2)+
      geom_line(aes(y=sM), color=sc[1], size=0.8)+
      scale_x_continuous(limits=c(0, MNE_max),
                         expand=c(0.0175,0.0175))+ 
      scale_y_continuous(limits=ylims, expand=c(0,0))+
      labs(x="TMNE", y="Exaggeration ratio")+
      annotate("text",label=plotlets[which(plotRels==x)], 
               x=MNE_max, y=ylims[2] - 0.02*(ylims[2]-ylims[1]), 
               hjust=1, vjust=1, size=6)+
      annotate("text",label=truevalString, x=xtext, 
               y=ylims[2] - 0.03*(ylims[2]-ylims[1]), 
               vjust=1, hjust=0, size=3.5, parse=TRUE)+
      annotate("text",label=bones2, x=MNE_max*0.65, 
               y=ylims[2] - 0.03*(ylims[2]-ylims[1]),
               vjust=1, size=3.5)+
      theme(panel.background=element_blank(),
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(),
            axis.title.x=eval(xtitle),
            axis.title.y=eval(ytitle),
            panel.border=element_rect(color="black", 
                                      fill=NA, size=0.6))

    return(explot)
    
  })
  
  # Create plot legend
  us <- runif(10, 7, 7.5)
  ls <- runif(10, 2, 2.5)
  ms <- runif(10, 4.9, 5.1)
  xs <- seq(0.4, 2.4, length.out=10)
  kdf <- data.frame(us=us, ls=ls, ms=ms, x1=xs, x2=xs+4.7)
  
  plotkey <- ggplot(data=kdf)+
    geom_ribbon(aes(ymin=ls, ymax=us, x=x1), fill=sc[2], alpha=0.2)+
    geom_line(aes(y=ms, x=x1), color=sc[2], size=0.8)+
    geom_ribbon(aes(ymin=ls, ymax=us, x=x2), fill=sc[1], alpha=0.2)+
    geom_line(aes(y=ms, x=x2), color=sc[1], size=0.8)+
    annotate("text", label="All correlations", color=sc[2], size=5,
             hjust=0, x=2.5, y=5)+
    annotate("text", label="p < 0.05 correlations", 
             color=sc[1], size=5, hjust=0, x=7.2, y=5)+
    scale_x_continuous(limits=c(0,10),expand=c(0,0))+
    scale_y_continuous(limits=c(0,10),expand=c(0,0))+
    theme_void()
  
  explots[[length(explots) + 1]] <- plotkey
  
  return(explots)
  
}


SignPlot <- function(sDF, bones, truevalDF, plotlets, 
                     xlabbool, ylabbool, max_y){
  
  # Remove simulations for which sign errors
  # were not calculated due to high MNE values or 
  # "None" relationship simulations
  sDF <- sDF[which(!is.na(sDF$r_med)), ]
  cors <- colnames(truevalDF)[3:5]
  
  relats <- truevalDF$str
  max_MNE <- max(sDF$MNE)
  
  # Create list of data frames for plotting
  pDFlist <- lapply(cors, function(z){
    
    nsigstat <- sDF[,which(colnames(sDF)==paste0("sign_", z))]
    sigstat <- sDF[,which(colnames(sDF)==paste0("sign_sig_", z))]
    df <- data.frame(MNE=rep(sDF$MNE, 2),
                     Relationship=rep(sDF$Relationship, 2),
                     sig=c(rep("n", nrow(sDF)), 
                           rep("y", nrow(sDF))),
                     y=c(nsigstat, sigstat),
                     stringsAsFactors=FALSE)
    df$g <- paste0(df$Relationship, df$sig)
    return(df)
  })

  # Create strings for plot text
  corlabs <- sapply(3:5, function(z){
    truevalString <- as.character(round(truevalDF[,z], 3))
    for(q in 1:length(truevalString)){
      if(nchar(truevalString[q] < 5)){
        zeros <- paste0(rep("0", (5-nchar(truevalString[q]))), 
                        collapse="")
        truevalString[q] <- paste0(truevalString[q], zeros)
      }
      
      tempname <- colnames(truevalDF[z])
      if(tempname=="beta"){
        truevalString[q] <- paste0("italic(beta)=='", 
                                  truevalString[q], "'")
      }else{
        truevalString[q] <- paste0("italic('", tempname,"')=='",
                                  truevalString[q], "'")
      }

    }
    return(truevalString)
  })
  
  # Positions and colors for plotting geometry
  xpos <- c(max_MNE*0.63, max_MNE*0.65, max_MNE*0.75,
            max_MNE*0.77, max_MNE*0.87)
  ypos <- c(max_y*0.83, max_y*0.74, max_y*0.67, max_y*0.60)
  gcols <- c("gray91", "gray68", "gray46",
             "lightgoldenrod", "orange1", "darkorange3")
  gsizes <- c(2, 1.2, 0.4)
  
  signplots <- lapply(1:length(pDFlist), function(z){
    
    xtitle <- ifelse(xlabbool[z], 
                     parse(text="element_text()"),
                     parse(text="element_blank()"))
    ytitle <- ifelse(ylabbool[z], 
                     parse(text="element_text()"),
                     parse(text="element_blank()"))
    
    if(bones=="26 parts"){
      bones2 <- "26 part types"
    }else{
      bones2 <- "6 part types"
    }
    
    signplot <- ggplot(pDFlist[[z]], aes(x=MNE, y=y))+
      geom_line(aes(group=g, color=g, size=Relationship,
                    linetype=sig))+
      annotate("text", label=plotlets[z], size=4.5,
               x=max_MNE, y=max_y*0.98, hjust=1, vjust=1)+
      annotate("text", label=bones2, size=4,
               x=xpos[3], y=max_y*0.96, vjust=1)+
      annotate("text", label=c("All correlations", 
                               "p < 0.05"),
               size=4, color=gcols[c(2, 5)],
               x=xpos[3:4], y=ypos[1],
               hjust=c(1, 0))+
      annotate("text", label=corlabs[,z], size=3.5,
               x=rep(xpos[1], 3), y=ypos[2:4],
               hjust=rep(1, 3), parse=TRUE)+
      annotate("segment", x=c(rep(xpos[2], 3), 
                              rep(xpos[4], 3)),
               xend=c(rep(xpos[3], 3), rep(xpos[5], 3)),
               y=rep(ypos[2:4], 2), yend=rep(ypos[2:4], 2),
               size=rep(gsizes, 2), color=gcols,
               linetype=c(rep("2111", 3), rep("solid", 3)))+
      scale_color_manual(values=c("Strongn"=gcols[1], 
                                "Moderaten"=gcols[2],
                                "Weakn"=gcols[3],
                                "Strongy"=gcols[4],
                                "Moderatey"=gcols[5],
                                "Weaky"=gcols[6]))+
      scale_size_manual(values=c("Strong"=gsizes[1],
                               "Moderate"=gsizes[2],
                               "Weak"=gsizes[3]))+
      scale_linetype_manual(values=c("y"="solid",
                                     "n"="2111"))+
      labs(x="TMNE", y="Sign error rate")+
      scale_x_continuous(limits=c(0, max_MNE),
                         expand=c(0.0175,0.0175))+ 
      scale_y_continuous(limits=c(0, max_y), expand=c(0,0))+
      theme(panel.background=element_blank(),
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            legend.position="none",
            axis.title.x=eval(xtitle),
            axis.title.y=eval(ytitle),
            panel.border=element_rect(color="black", 
                                      fill=NA, size=0.6))
 
    return(signplot)
    
  })
  
  return(signplots)
  
}


# Function to plot sampled Bayesian distributions
BayesPlot <- function(bayesDF, bayesDF2, cDF, bones, x_max, n_dists,  
                      relats, xlabbool, ylabbool, Plotlets){
  
  # Set y axis limits
  ybreaks <- seq(-0.036, 0.036, by=0.018)
  
  # Get positions of simulated MNE values closest
  # to each position in the vector of MNE values for
  # which to plot distributions
  maxwidth <- x_max/n_dists
  distbreaks <- seq(maxwidth/2, x_max-(maxwidth/2), 
                    length.out=n_dists)
  distMNE <- sapply(distbreaks, function(x){
    near <- abs(bayesDF$MNE - x)
    return(bayesDF$MNE[which.min(near)][1])
  })
  
  bayesDF <- bayesDF[which(bayesDF$MNE %in% distMNE), ]
  
  # Return maximum density value between the 26 parts
  # and 6 long bones data frames
  maxdens <- max(c(bayesDF$bu_den, bayesDF2$bu_dens))
  
  # Format density values in terms of MNE for plotting.
  denstandardized <- maxwidth*(bayesDF$bu_den/maxdens)
  bayesDF$densmax <- bayesDF$MNE + 0.5*denstandardized
  bayesDF$densmin <- bayesDF$MNE - 0.5*denstandardized
  
  # Create grouping variable
  bayesDF$g <- ifelse(bayesDF$dist=="prior",
                      paste0("a", bayesDF$MNE),
                      paste0("b", bayesDF$MNE))
  
  bayesplots <- lapply(1:length(relats), function(x){
    
    pbDF <- bayesDF[which(bayesDF$Relationship==relats[x]), ]
    
    # Get true value for relationship and format it for
    # plotting display
    if(relats[x]!="None"){
      truevalue <- cDF$beta[which(cDF$str==relats[x])]
      truevalString <- as.character(round(truevalue, 3))
      truelinealpha <- 1
      if(nchar(truevalString < 5)){
        zeros <- paste0(rep("0", (5-nchar(truevalString))), collapse="")
        truevalString <- paste0(truevalString, zeros)
      }
      truevalString <- paste0("beta=='", truevalString, "'")
    }else{
      truevalue <- truelinealpha <- 0
      truevalString <- paste0("beta=='0.000'")
    }
    
    # Conditional axis title display
    xtitle <- ifelse(xlabbool[x], 
                     parse(text="element_text()"),
                     parse(text="element_blank()"))
    ytitle <- ifelse(ylabbool[x], 
                     parse(text="element_text()"),
                     parse(text="element_blank()"))
    
    if(bones=="26 parts"){
      bones2 <- "26 part types"
    }else{
      bones2 <- "6 part types"
    }
    
    # Create plot
    bayesplot <- ggplot(pbDF, aes(xmin=densmin, xmax=densmax))+
      annotate("segment",x=-Inf, xend=Inf, y=0, yend=0, size=0.5)+
      geom_ribbon(aes(group=g, y=bu_val, fill=dist, alpha=dist))+
      annotate("segment",x=-Inf, xend=Inf, y=truevalue, yend=truevalue, 
               size=0.7, linetype=2, alpha=truelinealpha)+
      scale_fill_manual(values=c("prior"="grey", "post"="orchid3"))+
      scale_alpha_manual(values=c("prior"=0, "post"=1))+
      scale_x_continuous(limits=c(0, x_max), expand=c(0.0175,0.0175))+ 
      scale_y_continuous(limits=1.1*c(min(ybreaks), max(ybreaks)),
                         breaks=ybreaks, expand=c(0,0))+
      labs(x="TMNE", y=expression(beta))+
      annotate("text",label=Plotlets[x], 
               x=x_max, y=min(ybreaks), hjust=1, vjust=0, size=6)+
      annotate("text",label=truevalString, x=x_max*0.65, 
               y=min(ybreaks)+0.0035*(max(ybreaks)-min(ybreaks)), 
               vjust=0, hjust=0, size=4, parse=TRUE)+
      annotate("text",label=bones2, x=x_max*0.6, 
               y=min(ybreaks)+0.005*(max(ybreaks)-min(ybreaks)),
               vjust=0, hjust=1, size=4)+
      theme(panel.background=element_blank(),
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(), 
            axis.title.y=eval(xtitle),
            axis.title.x=eval(ytitle),
            panel.border=element_rect(color="black", 
                                      fill=NA, size=0.6),
            legend.position="none")
    
    return(bayesplot)
    
  })
  
  return(bayesplots)
  
}
