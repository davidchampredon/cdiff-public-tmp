###
###  RUN MONTE CARLO SIMULATIONS
###

R.library.dir    <- '../Rlibrary/lib'
library(dplyr);library(tidyr) #library(tidyverse)
library(snowfall)
library(abmcdiff,lib.loc = R.library.dir)
#library(profvis)

source('utils.R')
source('rooms-builder.R')

tictac <- as.numeric(Sys.time())

# Argument from the commande line:
args <- commandArgs(trailingOnly = TRUE)


#' DELETE? ALREADY IN "utils.R"
#' One monte carlo simulation.
# run_one <- function(seed_i, hcs, sim.prm, disease.prm, mvt_hcs) {
#     # IMPORTANT : do not forget "set.seed()" in this R script    
#     sim.prm   <- c(sim.prm, list(seedRNG  = seed_i)) 
#     one.simul <- abmcdiff_one_simul(hcs, sim.prm, disease.prm, mvt_hcs)
#     return(one.simul)
# }


# ---- RUN ----

run.comp <- TRUE

if(run.comp){
    
    
    # args <- c("prm_vax_strategy_0.csv","prm_vax_p_0.csv")
    if(length(args)<2){
        message('NOT ENOUGH ARGUMENTS!')
        stop()
    }
    
    fname.vax.strategy <- args[1] # fname.vax.strategy = "prm_vax_strategy_0.csv"
    fname.vax          <- args[2] # fname.vax = "prm_vax_p_0.csv"
    msg <- paste('>>> Running with:', 
                 fname.vax.strategy,
                 fname.vax,
                 '...')
    message(msg); print(msg)
    tmp <- run_simulation_mc_new(file.vax.strategy = fname.vax.strategy,
                             file.vax = fname.vax)
    
    
    #DEBUG:
    file.room.bed = 'prm_room_bed.csv'
    file.hosp.epi = 'prm_hosp.csv'
    file.vax      = 'prm_vax.csv'
    file.prm.simul = 'prm_sim.csv'
    file.disease = 'prm_disease.csv'
    rngseed = 123456
}
