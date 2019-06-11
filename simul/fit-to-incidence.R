###
###  RUN MONTE CARLO SIMULATIONS
###

R.library.dir    <- '../Rlibrary/lib'
library(tidyverse)
library(snowfall)
library(abmcdiff,lib.loc = R.library.dir)
library(profvis)

source('utils.R')
source('rooms-builder.R')

set.seed(1234)

tictac <- as.numeric(Sys.time())


do.run <- TRUE

if(do.run){
    
    # ENTER TARGET INCIDENCE HERE:
    target.inc.pd.10k <- 5.0
    
    # Define ABC priors:
    
    n.abc <- 3
    single.prm.fit <- TRUE
    
    if(single.prm.fit){
        # Single parameter fit:
        prm.to.fit <- c('sIdx_newPatient_prm1')
        abc.smp <- runif(n = n.abc, 
                         min = 0.005, 
                         max = 0.06)
        assign(x=prm.to.fit[1], value = abc.smp)
        message('Single parameter fit')
        message('ABC sampled values (first 50): ')
        message(paste(abc.smp[1:50], collapse =' ; '))
    }
    
    if(!single.prm.fit){ 
        # For multiple parameters fit
        prm.to.fit <- c('prev_new_admission',
                        'visit_room_halftime',
                        'contact_indiv_halftime')
        assign(x=prm.to.fit[1], value = runif(n = n.abc, min = 0.01, max = 0.09))
        assign(x=prm.to.fit[2], value = runif(n = n.abc, min = 0.01, max = 10))
        assign(x=prm.to.fit[3], value = runif(n = n.abc, min = 0.01, max = 10))
    }
    
    # Base parameters:
    
    param.base <- read.csv(file = 'prm_hosp.csv', 
                           stringsAsFactors = F, strip.white = TRUE)
    
    fname.fit <- 'prm_hosp_FIT.csv'
    
    mean_inc_pd <- function(df.ipd) {
        yrmax <- max(df.ipd$year)
        mean.incpd <- mean(df.ipd$inc.per.10000.pd[df.ipd$year==yrmax])
        return(mean.incpd)
    }
    
    incpd.vec <- vector(length = n.abc)
    
    for(i in 1:n.abc) {# 1:n.abc
        
        # Overwrite with ABC values:
        param.i <- param.base
        for(k in 1:length(prm.to.fit)){
            param.i$value[param.base$name == prm.to.fit[k]] <- as.numeric(get(prm.to.fit[k])[i])
        }
        
        # Create temporary file
        write.csv(x = param.i, 
                  file = fname.fit, 
                  quote = FALSE, row.names = FALSE)
        
        # Run simulation:
        message(paste('Fit ABC iter',i,'/',n.abc,'...'))
        tmp <- run_simulation_mc_new(file.vax.strategy = 'prm_vax_strategy_0.csv',
                                     file.vax          = 'prm_vax_p_0.csv',
                                     file.prm.simul    =  'prm_sim_fit.csv',
                                     file.hosp.epi     = fname.fit)
        
        # Calculate mean incidence per patient days:
        df <- tmp[['df']]
        df.ipd <- incidence_patientdays(df)
        incpd.vec[i] <- mean_inc_pd(df.ipd)
    }
    
    # Create dataframe of distance and prm values:
    df.dist <- data.frame(iter = 1:n.abc,
                          incpd = incpd.vec)
    for(j in 1:length(prm.to.fit)){
        x <- get(prm.to.fit[j])
        df.dist <- cbind(df.dist, x)
        names(df.dist)[ncol(df.dist)] <- prm.to.fit[j]
    }
    
    # Distance from target incidence:
    df.dist$dist <- abs(df.dist$incpd - target.inc.pd.10k)
    head(df.dist)
    tail(df.dist)
    
    # Identify the best prm set:
    idx.best <- which.min(df.dist$dist)
    prm.best <- df.dist[idx.best,]
    
    # Plots
    df.plot <- df.dist[,-c(1:2)]
    pairs(df.plot)
    
    # Display and save
    print(t(prm.best))
    save(list = c('df.dist', 'prm.best'),
         file = 'fit-to-incidence.RData')
    
}
