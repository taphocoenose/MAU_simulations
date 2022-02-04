### Skeletal parts and utility indices ###

# Code authored by Ryan Breslawski (rbreslawski@smu.edu)
# Last edited in R v 4.0.2 on a Windows 10 machine, Mar 2, 2021

# These functions are required for and called within SEAexamples.R. 
# They rely on functions from the rethinking, ggplot2, and patchwork 
# packages, which are called in SEAexamples.R.


# Function to fit a Bayesian Poisson regression model to 
# an input data frame and extract and summarize the model
# parameters
BayesModel <- function(df, MNE_name, Predictor_name, AF_name){
  
  if(length(Predictor_name)==1){
    
    modelname <- paste0(MNE_name, Predictor_name, collapse="")
    
    # Create data frame columns for model based on user input
    # arguments
    d <- list(MNE=df[ , which(colnames(df)==MNE_name)],
              U=df[ , which(colnames(df)==Predictor_name)],
              os=log(df[ , which(colnames(df)==AF_name)]),
              Ind=1:nrow(df))
    
    # Fit Bayesian model to sample
    m<- map(alist(
      MNE ~ dpois(lambda),
      log(lambda) <- os + a + bu*U,
      a ~ dnorm(0, 50),
      bu ~ dnorm(0, 0.01)), data=d,
      start=list(bu=0, a=0))
    
    # Extract samples and generate density values
    samples <- extract.samples(m, n=1e5)
    densn <- 500
    samplesdens1 <- density(samples$bu, n=densn)
    samplesdf <- data.frame(densID=rep(Predictor_name, densn),
                            x=samplesdens1$x, 
                            y=samplesdens1$y,
                            stringsAsFactors=FALSE)
    
    # Create a data frame for plotting model predictions
    simvals <- seq(0, 100, by=0.5)
    simdf <- t(sapply(simvals, function(x){
      
      medpred <- exp(samples$a + samples$bu*x)
      medpredsum <- quantile(medpred, probs=c(0.5, 0.025, 0.975))
      
      simsamples <- sapply(1:length(medpred), function(y) rpois(1, medpred[y]))
      
      samplessum <- quantile(simsamples, probs=c(0.025, 0.975))
      
      return(c(x, medpredsum, samplessum))
    }))
    
    simdf <- as.data.frame(simdf)
    colnames(simdf) <- c("x", "m_me", "m_lo", "m_up",
                         "pre_lo", "pre_up")
    
    fitsum <- data.frame(par=c("a", "bu"), med=rep(NA, 2), 
                         l025=rep(NA, 2), u975=rep(NA, 2), 
                         stringsAsFactors=FALSE)
    fitsum[1, 2:4] <- quantile(samples$a, c(0.5, 0.025, 0.975))
    fitsum[2, 2:4] <- quantile(samples$bu, c(0.5, 0.025, 0.975))
    
  }else{
    
    modelname <- paste0(MNE_name, Predictor_name[1], 
                        Predictor_name[2], collapse="")
    
    # Create data frame columns for model based on user input
    # arguments
    d <- list(MNE=df[ , which(colnames(df)==MNE_name)],
              U1=df[ , which(colnames(df)==Predictor_name[1])],
              U2=df[ , which(colnames(df)==Predictor_name[2])],
              os=log(df[ , which(colnames(df)==AF_name)]),
              Ind=1:nrow(df))
    
    # Fit Bayesian model to sample
    m <- map(alist(
      MNE ~ dpois(lambda),
      log(lambda) <- os + a + bu1*U1 + bu2*U2,
      a ~ dnorm(0, 50),
      bu1 ~ dnorm(0, 0.01),
      bu2 ~ dnorm(0, 0.01)), data=d,
      start=list(bu1=0, bu2=0, a=0))
    
    # Extract samples and generate density values
    samples <- extract.samples(m, n=1e5)
    densn <- 500
    samplesdens1 <- density(samples$bu1, n=densn)
    samplesdens2 <- density(samples$bu2, n=densn)
    samplesdf <- data.frame(densID=c(rep(Predictor_name[1], 
                                         densn),
                                     rep(Predictor_name[2], 
                                         densn)),
                            x=c(samplesdens1$x, samplesdens2$x),
                            y=c(samplesdens1$y, samplesdens2$y),
                            stringsAsFactors=FALSE)
    
    simdf <- NA
    
    fitsum <- data.frame(par=c("a", "bu1", "bu2"),
                         med=rep(NA, 3), l025=rep(NA, 3),
                         u975=rep(NA, 3), stringsAsFactors=FALSE)
    fitsum[1, 2:4] <- quantile(samples$a, c(0.5, 0.025, 0.975))
    fitsum[2, 2:4] <- quantile(samples$bu1, c(0.5, 0.025, 0.975))
    fitsum[3, 2:4] <- quantile(samples$bu2, c(0.5, 0.025, 0.975))
  }
  
  return(list(MNE_vals=MNE_name, predictor=Predictor_name,
              data=d, samplesdf=samplesdf, summary=fitsum, 
              plotdata=simdf, model=m))
}


# Function to plot the linear relationship with posterior
# predictions for fitted models with a single predictor
PostModelPlot <- function(x, plotlet){
  
  d <- x$plotdata
  d$U <- d$x
  x$data$MAU <- x$data$MNE/exp(x$data$os)
  
  y_max <- max(c(d$m_up, x$data$MAU))
  
  # Truncate the predictive distribution if it extends
  # above the plot y range
  d$pre_up[which(d$pre_up > y_max)] <- y_max
  
  p <- ggplot(d, aes(x=U))+
    geom_ribbon(aes(ymin=pre_lo, ymax=pre_up), fill="gray81")+
    geom_line(aes(y=m_me), color="white", size=0.8)+
    geom_line(aes(y=m_lo), color="white", linetype="dashed")+
    geom_line(aes(y=m_up), color="white", linetype="dashed")+
    annotate("point", y=x$data$MAU, x=x$data$U)+
    annotate("text", label=plotlet, x=100, y=0,
             hjust=1, vjust=0)+
    labs(x=x$predictor, y="MAU")+
    scale_x_continuous(limits=c(0, 100), expand=c(0.0175, 
                                                  0.0175))+ 
    scale_y_continuous(limits=c(0, y_max), expand=c(0.02,0.1))+
    theme(panel.background=element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.border=element_rect(color="black", 
                                    fill=NA, size=0.6))
  
    return(p)
}


# Function to plot the posterior parameter distributions
# across models
PostParsPlot <- function(x, modcols, plotlet, xtitlebool, 
                         ltitle, ldisply, mnames=NA){
  
  # Combine data frames across models
  d <- lapply(1:length(x), function(y){
    dx <- x[[y]]$samplesdf
    dx$Model <- paste(x[[y]]$predictor, collapse="+")
    return(dx)
  })
  d <- do.call("rbind", d)
  
  # If models are distinguished by datasets rather than
  # predictors, rename models in the d data drame
  if(length(unique(d$Model))==1){
    sn <- nrow(x[[1]]$samplesdf)
    d$Model <- c(sapply(mnames, function(y) rep(y, sn)))
  }
  
  d$g <- paste0(d$Model, d$densID)
  
  # Get parameter names
  pars <- unique(d$densID)
  
  # Model colors for plotting
  models <- unique(d$Model)
  names(modcols) <- models
  
  # Rescale and reposition densities for plotting
  d$y <- d$y/max(d$y)
  d$yy <- sapply(d$densID, function(y) which(pars==y))
  d$y_min <- d$yy - d$y/2
  d$y_max <- d$yy + d$y/2
  
  # Remove densities that are so small they are not
  # well observed on the plot
  d <- d[which(d$yy > 0.001), ]
  
  x_min <- min(d$x) - 0.2*(max(d$x) - min(d$x))
  
  xtitle <- ifelse(xtitlebool, 
                   parse(text="element_text()"),
                   parse(text="element_blank()"))
  
  # Reformat par names for plotting
  pars <- gsub("_", " ", pars)
  
  p <- ggplot(d, aes(x=x, ymin=y_min, ymax=y_max))+
    geom_ribbon(aes(group=g, fill=Model), alpha=0.6)+
    scale_fill_manual(values=modcols)+
    annotate("text", label=pars, x=x_min,
             y=1:length(pars), hjust=0)+
    annotate("text", label=plotlet, y=0.5, x=max(d$x), 
             vjust=0, hjust=1)+
    labs(x=expression(beta), fill=ltitle)+
    theme(panel.background=element_blank(),
          panel.grid.major.y=element_blank(), 
          panel.grid.major.x=element_line(color="gray89"),
          panel.grid.minor=element_blank(),
          legend.position=ldisply,
          axis.title.y=element_blank(),
          axis.title.x=eval(xtitle),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.border=element_rect(color="black", 
                                    fill=NA, size=0.6))
  
  return(p)
  
}
