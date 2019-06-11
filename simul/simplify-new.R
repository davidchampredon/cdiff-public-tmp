###
###   SIMPLIFY THE DATAFRAME OUTPUT FROM SIMULATIONS
###   (will be typically used for downstream economic analysis)
###

library(dplyr);library(tidyr) #library(tidyverse)
source('utils.R')

# ---- Read RData files ----

rdatafiles <- system('ls simul-mc*.RData', intern = TRUE)
rdatafiles <- rdatafiles[!grepl('simpl',rdatafiles)]

if(length(rdatafiles)==0){
    message('ERROR: RData files not found. Aborting.')
    stop()
}

# ---- Read Vaccine strategy:
vaxstrategy <- read.csv('prm_vax_strategy_frailty.csv', as.is = TRUE)
vax.start.time <- as.numeric(vaxstrategy$value[vaxstrategy$name == "vax_start_time"])

# ---- Functions ----

#' Simplify the output dataframe.
simplify_df <- function(df, vax.start.time) {
    
    # NOTE:
    # Not using `tidyverse` functions (filter(), mutate(),...)
    # because it's too slow for those huge dataframes:
    
    # Interested in CDI cases or vaccinated patients only:
    # tmp <- df[df$symptomatic | df$isVax,]
    tmp <- df[df$wasColonized | df$isVax,]
    # tmp = df
    
    # Only interested in what happens after vaccination:
    tmp <- tmp[tmp$time > vax.start.time,]

    # Makes a key to untangle UIDs and MC:
    tmp$key <- paste(tmp$mc, tmp$uid, sep='_')
    
    # Multiple uids, keep the latest (it has the most information):
    tmp = tmp %>%
        group_by(key) %>%
        top_n(1,time)
        
    res = tmp[,c('mc', 
                 'uid', 
                 'timeAdmission',
                 'totDurSymptoms',
                 'isoDurFirst', 
                 'isoDurRelapses',
                 'isVax', 
                 'wasColonized',
                 'colectomy') ]
    
    return(res)
}


check_simplified <- function(a, fn='') {   #DEBUG:  a = simple.df    a=res
    print('\n  --- Checking simplified dataframe :')
    print(fn)
    
    # rows of symptomatic cases only"
    idx.sympt = a$totDurSymptoms>0
    
    print(paste('Total MC iter: ',length(unique(a$mc))))
    print(paste('Proportion vax among cases:',mean(a$isVax[idx.sympt])))
    print(paste('Proportion collectomy:',mean(a$colectomy)))
    print('First isolation duration:')
    print(summary(a$isoDurFirst[idx.sympt]))
    print('Relapses isolation duration:')
    print(summary(a$isoDurRelapses[a$isoDurRelapses>0]))
    
    print("Mismatch [symptoms - isolation] total durations:")    
    z = a$totDurSymptoms[idx.sympt] - (a$isoDurFirst[idx.sympt] + a$isoDurRelapses[idx.sympt])
    print(summary(z))
    
    n.vax.col <- nrow(a[a$isVax==TRUE & a$wasColonized==TRUE, ])
    n.sympt   <- sum(idx.sympt)
    prop.sympt.vax.col =  round(n.sympt/n.vax.col,4)
    print(paste("Proportion symptomatic among vaccinated and colonized:", 
                prop.sympt.vax.col))
    
    print('  --- Checks completed ---')
}

simple_rdatafile_name <- function(x){
    y <- strsplit(x,split = '.', fixed = TRUE)[[1]]
    z <- paste0(y[1],'-simplified.RData')
    return(z)
}


# ---- Run ----

for(x in rdatafiles){  #DEBUG: x=rdatafiles[1]
    message(paste('Simplifying',x,'...'))
    load(x)
    message('Rdata file loaded.\n Symplifying...')
    simple.df <- simplify_df(df, vax.start.time)

    check_simplified(simple.df,x)
    message('\n Symplification done. Saving simple RData file...')
    save(list = 'simple.df', 
         file = simple_rdatafile_name(x))
    rm(df)
    message('DONE.')
}

message('\n  ---- All RData files simplified ---- \n')
