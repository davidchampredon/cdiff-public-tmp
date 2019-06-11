###
###  Calculate sensitivity to all parameters
###

R.library.dir    <- '../Rlibrary/lib'
library(tidyverse)
library(snowfall)
library(abmcdiff,lib.loc = R.library.dir)

source('utils.R')
source('rooms-builder.R')

set.seed(1234)

tictac <- as.numeric(Sys.time())

do.run <- TRUE

# ==== Functions ====

# Calculate mean incidence per patient days
mean_inc_pd <- function(df.ipd) {
    yrmax <- max(df.ipd$year)
    mean.incpd <- mean(df.ipd$inc.per.10000.pd[df.ipd$year==yrmax])
    return(mean.incpd)
}

# Retrieve response variable from simulations:
get_resp_var <- function(sim) {
    return(sim[['df']] %>%
               incidence_patientdays() %>% 
               mean_inc_pd()
    )
}

# ==== RUN ====

if(do.run){
    # Base parameters:
    filename.prm <- 'prm_hosp.csv'
    param.base <- read.csv(file = filename.prm, 
                           stringsAsFactors = F, 
                           strip.white = TRUE)
    fname.tmp <- paste0('tmp_',filename.prm)
    
    # Filter out non numerical values:
    idx.keep <- !is.na(as.numeric(param.base$value))
    param.string <- param.base[!idx.keep,]
    param.base <- param.base[idx.keep,]
    param.base$value <- as.numeric(param.base$value)
    param.base
    param.bmp <- param.base
    
    # Reference simulation: 
    
    message(paste('Calculating reference value for sensitivities...'))
    simref <- run_simulation_mc_new(file.vax.strategy = 'prm_vax_strategy_0.csv',
                                 file.vax          = 'prm_vax_p_0.csv',
                                 file.prm.simul    = 'prm_sim_sensi.csv',
                                 file.hosp.epi     = filename.prm)
    
    # Reference response variable: (incidence per patient days)
    ref.value <- get_resp_var(simref)
    
    # Bumped parameter values
    bmp <- 0.1
    N <- nrow(param.base)
    sensi <- rep(NA,N)
    val.bmp <- rep(NA,N)
    
    # ----------------------- DEBUG
    # N = 2
    # ----------------------- DEBUG
    
    for(i in 1:N){    # i=1
        
        # Overwrite with bumped values:
        if(param.base$value[i] > 1e-2) val.bmp[i]  <- param.base$value[i] * (1+bmp)
        if(param.base$value[i] <= 1e-2) val.bmp[i]  <- param.base$value[i] * 10
        if(param.base$value[i] <= 1e-6) val.bmp[i]  <- 1e-3
        param.bmp$value[i] <- val.bmp[i]
        
        # Create temporary file
        write.csv(x = rbind(param.string,param.bmp), 
                  file = fname.tmp, 
                  quote = FALSE, 
                  row.names = FALSE)
        
        # Run simulation:
        message(paste('Calculating sensitivity',i,'/',N,':',param.base$name[i]))
        tmp <- run_simulation_mc_new(file.vax.strategy = 'prm_vax_strategy_0.csv',
                                     file.vax          = 'prm_vax_p_0.csv',
                                     file.prm.simul    = 'prm_sim.csv',
                                     file.hosp.epi     = fname.tmp)
        
        # Calculate sensitivity:
        sensi[i] <- get_resp_var(tmp) - ref.value
        
        # Restore initial value:
        param.bmp$value[i] <- param.base$value[i]
    }
    
    df.sensi <- data.frame(name = param.base$name, 
                           sensi = sensi,
                           init_value = param.base$value, 
                           bumped_value = val.bmp)
    
    save.image(file = 'sensi.RData')
    
    df.sensi$name.o <- factor(df.sensi$name, 
                              levels = df.sensi$name[order(df.sensi$sensi)])
    g <- df.sensi %>%
        ggplot()+
        geom_point(aes(x=name.o, y=sensi),
                   size=2)+
        coord_flip()
    plot(g)
    
}
