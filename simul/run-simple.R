###
###  RUN MONTE CARLO SIMULATIONS
###

R.library.dir    <- '../Rlibrary/lib'
library(tidyr)
library(dplyr)
library(ggplot2)
library(snowfall)
library(abmcdiff,lib.loc = R.library.dir)

source('utils.R')
source('rooms-builder.R')

tictac <- as.numeric(Sys.time())


do.run <- TRUE

if(do.run){
    
    fname <- c('prm_vax_strategy_0.csv')
    # fname <- c('prm_vax_strategy_frailty.csv')
    
    for(i in 1:length(fname)){
        msg <- paste('>>> Running with:', fname[i],'...', i,'/',length(fname))
        message(msg); print(msg)
        tmp <- run_simulation_mc_new(file.vax.strategy = fname[i],
                                 file.room.bed = 'prm_room_bed.csv',
                                 file.hosp.epi = 'prm_hosp.csv',
                                 file.vax      = 'prm_vax_p_0.csv', #'prm_vax_p_0.csv',
                                 file.prm.simul= 'prm_sim.csv',
                                 file.disease  = 'prm_disease.csv',
                                 rngseed = 123456)
    }
    
    # --- DEBUG:
    file.room.bed = 'prm_room_bed.csv'
    file.hosp.epi = 'prm_hosp.csv'
    file.vax      = 'prm_vax_p_0.csv'
    file.prm.simul= 'prm_sim.csv'
    file.disease  = 'prm_disease.csv'
    rngseed = 123456
    # ---
}
