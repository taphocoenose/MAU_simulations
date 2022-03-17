### Skeletal parts and utility indices ###

# Code authored by Ryan Breslawski (rbreslawski@smu.edu)
# Last edited in R v 4.0.4 on a Windows 10 machine, Jan 25, 2022

# These functions are required for and called within SEAsim.R. They rely
# on functions from the rethinking and parallel packages, which are 
# called in SEAsim.R.

# Function to obtain and plot population correlations
PopBones <- function(x, bones, plotlets){
  
  PopPlots <- list()
  
  cDF <- data.frame(bones=rep(bones, 3), 
                    str=c("Strong", "Moderate", "Weak"),
                    r=rep(NA, 3), rho=rep(NA, 3), 
                    beta=rep(NA, 3))
  
  for(i in 4:6){
    
    x2 <- x
    x2$pMAU <- x2[,i]/x2$Anatomical.Frequency
    x2$pMAU <- 100*x2$pMAU/max(x2$pMAU)
    x2$pseudoCount <- floor(x2[,i]*1e6)
    
    r <- cor.test(x2$Utility, x2$pMAU, 
                  method="pearson")$estimate
    rho <- cor.test(x2$Utility, x2$pMAU, 
                    method="spearman", exact=F)$estimate
    beta <- coef(summary(glm(pseudoCount ~ Utility + offset(os),
                             data=x2, family=poisson())))[2,1]
    
    cDF[i-3, 3:5] <- c(r, rho, beta)
    
    plotrho <- as.character(round(rho, 3))
    if(nchar(plotrho)<5){
      zeros <- paste0(rep("0", 5-nchar(plotrho)), collapse="")
      plotrho <- paste0(plotrho, zeros)
    }
    plotr <- as.character(round(r, 3))
    if(nchar(plotr)<5){
      zeros <- paste0(rep("0", 5-nchar(plotr)), collapse="")
      plotr <- paste0(plotr, zeros)
    }
    
    xtitle <- ifelse(i==6, 
                     parse(text="element_text()"),
                     parse(text="element_blank()"))
    ytitle <- ifelse(bones=="26 parts", 
                     parse(text="element_text()"),
                     parse(text="element_blank()"))
    
    relname <- paste(cDF$str[i-3], "relationship")
    
    if(bones=="26 parts"){
      bones2 <- "26 part types"
    }else{
      bones2 <- "6 part types"
    }
    
    popplot <- ggplot(data=x2, aes(x=Utility, y=pMAU))+
      geom_smooth(method="lm", color="grey", size=1, 
                  fullrange=TRUE, se=FALSE)+
      geom_point(shape=16, color="black")+
      labs(x="MGUI",y="%MAU")+
      annotate("text", label=bones2, y=98, x=0, hjust=0, color="grey")+
      annotate("text",label=plotlets[i-3],x=max(x2$Utility),
               y=4,hjust=1,vjust=0,size=6)+
      annotate("text", label=relname, x=max(x2$Utility)*0.55, 
               y=17, hjust=0, color="grey")+
      annotate("text",label=paste("italic('rho')=='", plotrho, "'"), 
               x=max(x2$Utility)*0.55, y=5, hjust=0, parse=TRUE,
               color="grey")+
      scale_x_continuous(limits=c(0, 100), 
                         expand=c(0.0175,0.0175))+ 
      scale_y_continuous(limits=c(0,105), 
                         breaks=seq(from=0,to=100, by=25),
                         expand=c(0.0175,0.0175))+
      theme(panel.background=element_blank(),
            panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(),
            axis.title.x=eval(xtitle),
            axis.title.y=eval(ytitle),
            panel.border=element_rect(color="black", 
                                      fill=NA, size=0.6))
    
    PopPlots[[i-3]] <- popplot
    
  }
  
  return(list(PopPlots=PopPlots, cDF=cDF))
}

# Define function to resample MNE values and run correlations
# for a given parameter set
SimCors <- function(x3, nCor1, simMNE1, MSM, NPM, Bayes){
  
  if((simMNE1 > MSM) | (NPM > nCor1)){
    nCor1 <- NPM
  }
  
  # Sample elements for simulation
  boneSamples <- sample(x=x3$Element, size=nCor1*simMNE1, 
                        replace=T, prob=x3$prob)
  boneSamples <- matrix(boneSamples, nrow=nCor1)
  
  # Simulate for parameter set j
  simResultsJ <- lapply(1:nrow(boneSamples), function(z){
    
    MNE <- sapply(x3$Element, function(v){
      length(which(boneSamples[z, ]==v))
    })
    MAU <- MNE/x3$Anatomical.Frequency
    
    rcor <- cor.test(x3$Utility, MAU, method="pearson")
    if(is.na(rcor$estimate)){
      rcor$estimate <- 0
      rcor$p.value <- 1
    }
    
    rhocor <- cor.test(x3$Utility, MAU, method="spearman", exact=F)
    if(is.na(rhocor$estimate)){
      rhocor$estimate <- 0
      rhocor$p.value <- 1
    }
    
    # Data frame for GLMs
    d <- data.frame(MNE=MNE, os=x3$os, Utility=x3$Utility)
    
    beta <- glm(MNE ~ Utility + offset(os), data=d, family=poisson())
    beta <- coef(summary(beta))
    if(is.na(beta[2,1])){
      beta[2,1] <- 0
      beta[2,4] <- 1
    }
    
    rDF <- data.frame(r=rcor$estimate, r_p=rcor$p.value,
                      rho=rhocor$estimate, rho_p=rhocor$p.value,
                      beta=beta[2,1], beta_p=beta[2,4])
    
    if(Bayes){
      
      beta_Bayes <- map(alist(
        MNE ~ dpois(lambda),
        log(lambda) <- os + a + bu*Utility,
        a ~ dnorm(0,5),
        bu ~ dnorm(0,0.01)), data=d, 
        start=list(bu=0, a=0))@coef[1]
      
      if(is.na(beta_Bayes)){beta_Bayes <- 0}
      
      rDF$beta_Bayes[1] <- beta_Bayes
    }
    
    return(rDF)
  })
  
  returnedCors <- do.call("rbind", simResultsJ)
  # Convert p values to binary significance character
  for(w in c(2, 4, 6)){
    returnedCors[,w] <- ifelse(returnedCors[,w]>0.05, "n", "y") 
  }
  
  return(returnedCors)
}


# Define function for simulation over all parameter sets
SimBones <- function(x1, cDF, MNE_values, nCor, n_per_MNE, 
                     MNE_sig_max, rnames1){
  
  # Parameter sets over which to perform simulation
  simPars <- data.frame(Relationship=rep(rnames1, length(MNE_values)),
                        MNE=sort(rep(MNE_values, length(rnames1))),
                        stringsAsFactors=FALSE)
  
  # make cluster and export variables and rethinking
  cl <- makeCluster(detectCores()-1)
  clusterExport(cl, varlist=list("x1", "nCor", "SimCors", "cDF", 
                                 "simPars" ,"MNE_sig_max", "n_per_MNE"),
                envir=environment())
  clusterCall(cl, function() library(rethinking))
  
  simtext <- paste0("'", cDF$bones[1], "' simulation:")
  cat(paste("Begin", simtext, Sys.time(), "\n"))
  
  # Simulate over each parameter set
  simResults <- clusterApplyLB(cl, 1:nrow(simPars), function(j){
    
    relat <- simPars$Relationship[j]
    simMNE <- simPars$MNE[j]
    x2 <- x1
    colnames(x2)[which(colnames(x2)==relat)] <- "prob"
    
    # Simulate samples and correlations for parameter set J
    AllCors <- SimCors(x3=x2, nCor1=nCor, simMNE1=simMNE, 
                       MSM=MNE_sig_max, NPM=n_per_MNE,
                       Bayes=T)
    
    # Summarize simulation results for parameter set j
    sum_PropSig_r <- length(which(AllCors[,2]=="y"))/nCor
    sum_PropSig_rho <- length(which(AllCors[,4]=="y"))/nCor
    sum_PropSig_beta <- length(which(AllCors[,6]=="y"))/nCor
    
    # If not simulating for the "None" relationship and
    # MNE value is less than or equal to MNE_sig_max, 
    # get exaggeration ratio and sign error statistics
    if((relat != "None") & (simMNE <= MNE_sig_max)){
      
      # Known population correlations
      k_r <- cDF$r[which(cDF$str==relat)]
      k_rho <- cDF$rho[which(cDF$str==relat)]
      k_beta <- cDF$beta[which(cDF$str==relat)]
      
      # Exaggeration ratio stats
      pq <- c(0.5, 0.025, 0.975)
      sum_r <- quantile(abs(AllCors[,1])/k_r, probs=pq)
      sum_rho <- quantile(abs(AllCors[,3])/k_rho, probs=pq)
      sum_beta <- quantile(abs(AllCors[,5])/k_beta, probs=pq)
      sum_beta_Bayes <- quantile(abs(AllCors[,7])/k_beta, probs=pq)
      
      # Sign error stats
      sign_r <- length(which(AllCors[,1] < 0))/nCor
      sign_rho <- length(which(AllCors[,3] < 0))/nCor
      sign_beta <- length(which(AllCors[,5] < 0))/nCor
      sign_bayes <- length(which(AllCors[,7] < 0))/nCor
        
      # Store significant correlations from first set of simulations
      # in separate data frames
      SigCorsr <- AllCors[which(AllCors[,2]=="y"),1:2]
      SigCorsrho <- AllCors[which(AllCors[,4]=="y"),3:4]
      SigCorsbeta <- AllCors[which(AllCors[,6]=="y"),5:6]
      
      # Simulated correlations until nCors significant
      # correlations are reached
      NotReachedBool <- rep(T, 3)
      while(any(NotReachedBool)){
        
        SigCorsTrial <- SimCors(x3=x2, nCor1=nCor, simMNE1=simMNE, 
                                MSM=MNE_sig_max, NPM=n_per_MNE,
                                Bayes=F)
        
        if(NotReachedBool[1]){
          SigCorsrTrial <- SigCorsTrial[which(SigCorsTrial[,2]=="y"),1:2]
          SigCorsr <- rbind(SigCorsr, SigCorsrTrial)
          if(nrow(SigCorsr) >= nCor){
            NotReachedBool[1] <- F
            SigCorsr <- SigCorsr[1:nCor,]
          }
        }
        
        if(NotReachedBool[2]){
          SigCorsrhoTrial <- SigCorsTrial[which(SigCorsTrial[,4]=="y"),3:4]
          SigCorsrho <- rbind(SigCorsrho, SigCorsrhoTrial)
          if(nrow(SigCorsrho) >= nCor){
            NotReachedBool[2] <- F
            SigCorsrho <- SigCorsrho[1:nCor,]
          }
        }
        
        if(NotReachedBool[3]){
          SigCorsbetaTrial <- SigCorsTrial[which(SigCorsTrial[,6]=="y"),5:6]
          SigCorsbeta <- rbind(SigCorsbeta, SigCorsbetaTrial)
          if(nrow(SigCorsbeta) >= nCor){
            NotReachedBool[3] <- F
            SigCorsbeta <- SigCorsbeta[1:nCor,]
          }
        }
      }
      
      # Exaggeration ratio stats
      sum_sig_r <- quantile(abs(SigCorsr[,1])/k_r, probs=pq)
      sum_sig_rho <- quantile(abs(SigCorsrho[,1])/k_rho, probs=pq)
      sum_sig_beta <- quantile(abs(SigCorsbeta[,1])/k_beta, probs=pq)
      
      # Sign error stats
      sign_sig_r <- length(which(SigCorsr[,1] < 0))/nCor
      sign_sig_rho <- length(which(SigCorsrho[,1] < 0))/nCor
      sign_sig_beta <- length(which(SigCorsbeta[,1] < 0))/nCor
      
    }else{
      
      # If the simulation is for relationship "None" or for an
      # MNE value higher than MNE_sig_max, output NA values
      # for exaggeration ratio and sign error stats
      sum_r <- sum_rho <- sum_beta <- sum_beta_Bayes <- rep(NA, 3)
      sum_sig_r <- sum_sig_rho <- sum_sig_beta <- rep(NA, 3)
      sign_r <- sign_rho <- sign_beta <- sign_beta_Bayes <- NA
      sign_sig_r <- sign_sig_rho <- sign_sig_beta <- NA
      
    }
    
    sumDF <- data.frame(Relationship=relat, MNE=simMNE,
                        PropSig_r=sum_PropSig_r,
                        PropSig_rho=sum_PropSig_rho,
                        PropSig_beta=sum_PropSig_beta,
                        r_med=sum_r[1], 
                        r_lo=sum_r[2], 
                        r_up=sum_r[3],
                        rho_med=sum_rho[1], 
                        rho_lo=sum_rho[2], 
                        rho_up=sum_rho[3], 
                        beta_med=sum_beta[1], 
                        beta_lo=sum_beta[2], 
                        beta_up=sum_beta[3],
                        beta_Bayes_med=sum_beta_Bayes[1], 
                        beta_Bayes_lo=sum_beta_Bayes[2], 
                        beta_Bayes_up=sum_beta_Bayes[3],
                        sig_r_med=sum_sig_r[1], 
                        sig_r_lo=sum_sig_r[2], 
                        sig_r_up=sum_sig_r[3], 
                        sig_rho_med=sum_sig_rho[1], 
                        sig_rho_lo=sum_sig_rho[2], 
                        sig_rho_up=sum_sig_rho[3], 
                        sig_beta_med=sum_sig_beta[1], 
                        sig_beta_lo=sum_sig_beta[2], 
                        sig_beta_up=sum_sig_beta[3],
                        sign_r=sign_r,
                        sign_rho=sign_rho,
                        sign_beta=sign_beta,
                        sign_sig_r=sign_sig_r,
                        sign_sig_rho=sign_sig_rho,
                        sign_sig_beta=sign_sig_beta,
                        stringsAsFactors=FALSE)
    
    # Subsample all correlation samples DF for plotting. Also,
    # add fields for the simulation j parameters.
    if(nrow(AllCors) > n_per_MNE){
      AllCors <- AllCors[sample(1:nCor, n_per_MNE),]
    }
    AllCors$Relationship <- rep(relat, nrow(AllCors))
    AllCors$MNE <- rep(simPars$MNE[j], nrow(AllCors))
    
    
    # Extract one sample for parameter set J for plotting
    # a sample set of prior and posterior distributions
    boneSamples <- sample(x=x2$Element, size=simMNE, 
                          replace=T, prob=x2$prob)
    sMNE <- sapply(x2$Element, function(v){
      length(which(boneSamples==v))
    })
    sdata <- data.frame(MNE=sMNE, Utility=x2$Utility, os=x2$os)
    # Fit Bayesian model to sample
    smodel <- map(alist(
      MNE ~ dpois(lambda),
      log(lambda) <- os + a + bu*Utility,
      a ~ dnorm(0,5),
      bu ~ dnorm(0,0.01)), data=sdata, 
      start=list(bu=0, a=0))
    
    # Extract posterior and prior samples, generate densities,
    # and aggregate this density data into a data frame
    npts <- 1000
    bu_post <- density(extract.samples(smodel, n=1e4)$bu, n=npts)
    bu_prior <- density(rnorm(n=1e4, 0, 0.01), n=npts)
    sampleBayes <- data.frame(dist=c(rep("prior", npts), 
                                     rep("post", npts)),
                              MNE=rep(simMNE, 2*npts),
                              Relationship=rep(relat, 2*npts),
                              bu_val=c(bu_prior$x, bu_post$x),
                              bu_den=c(bu_prior$y, bu_post$y),
                              stringsAsFactors=FALSE)
    
    return(list(SimSummary=sumDF, sims=AllCors, sampleBayes=sampleBayes))
    
  })
  
  stopCluster(cl)
  
  # Extract 3 lists from SimResults. The first list contains
  # dataframes of simulation summaries. The second list contains
  # raw simulated correlations. The third list contains sample
  # prior and posterior distributions for plotting.
  SimsSummary <- lapply(1:length(simResults), function(o){
    simResults[[o]]$SimSummary
  })
  SimsSummary <- do.call("rbind", SimsSummary)
  SimsValues <- lapply(1:length(simResults), function(o){
    simResults[[o]]$sims
  })
  SimsValues <- do.call("rbind", SimsValues)
  SampleBayes <- lapply(1:length(simResults), function(o){
    simResults[[o]]$sampleBayes
  })
  SampleBayes <- do.call("rbind", SampleBayes)
  
  cat(paste("Completed", simtext, Sys.time(), "\n\n"))
  
  return(list(SimsSummary=SimsSummary, SimsValues=SimsValues,
              SampleBayes=SampleBayes))
  
}
