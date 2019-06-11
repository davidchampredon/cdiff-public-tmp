library(data.table)


create.list <- function(z) {
    param <- list()
    for(i in 1:nrow(z)){
        test.num <- as.numeric(z$value[i])
        print(paste(i,'test numeric: ',test.num))
        if(is.na(test.num)){
            param[[z$name[i]]] <- z$value[i]
            print(paste(i,'not numeric: ',z$name[i],z$value[i]))
        }
        else{
            param[[z$name[i]]] <- as.numeric(z$value[i])
            print(paste(i,'numeric: ',z$name[i],z$value[i]))
        }
    }
    return(param)
}

#' One monte carlo simulation.
run_one <- function(seed_i, hcs, sim.prm, disease.prm, mvt_hcs) {
    # IMPORTANT : do not forget "set.seed()" in this R script    
    
    do.debug <- FALSE
    if(do.debug){
        message(paste('DEBUG run_one()',seed_i))
        print(hcs)
        print(disease.prm)
    }
    
    sim.prm   <- c(sim.prm, list(seedRNG  = seed_i)) 
    one.simul <- abmcdiff_one_simul(hcs, sim.prm, disease.prm, mvt_hcs)
    return(one.simul)
}


run_simulation_mc <- function(file.vax.strategy,
                              file.room.bed = 'prm_room_bed.csv',
                              file.hosp.epi = 'prm_hosp.csv',
                              file.vax = 'prm_vax.csv',
                              file.prm.simul = 'prm_sim_fit.csv',
                              file.disease = 'prm_disease.csv',
                              rngseed = 123456) {
    
    set.seed(rngseed)
    
    # ---- Set up ----
    # --- Health Care Settings definitions:
    file.room.hosp <- 'rooms_hosp.csv'
    prm.room.bed <- create.list(read.csv(file = file.room.bed ,#'prm_room_bed.csv', 
                                         as.is = T,
                                         strip.white = T))
    
    r.hosp <- build_rooms(n.iso   = prm.room.bed$n.iso, 
                          n.wards = prm.room.bed$n.wards,
                          rooms.per.ward.mean = prm.room.bed$rooms.per.ward.mean, 
                          beds.per.room.mean  = prm.room.bed$beds.per.room.mean,
                          min.bed.per.room = prm.room.bed$min.bed.per.room,
                          max.bed.per.room = prm.room.bed$max.bed.per.room,
                          filename = file.room.hosp)
    summary_rooms(r.hosp)
    plot_rooms(r.hosp)
    
    prm.hosp <- full.param.HCS(file_infrastruct = file.room.hosp, 
                               file_epi         = file.hosp.epi) 
    
    # --- Vaccine features
    
    prm.vax <- create.list(read.csv(file = file.vax, #'prm_vax.csv', 
                                    as.is = T,
                                    strip.white = T))
    prm.vax.features <- list(vax_name = 'vaxificile',
                             vax_efficacy_lo = prm.vax$vax_efficacy_lo,
                             vax_efficacy_hi = prm.vax$vax_efficacy_hi,
                             vax_proba_prevent_colonization = prm.vax$vax_proba_prevent_colonization,
                             vax_halflife = prm.vax$vax_halflife)
    
    # --- Vaccination strategies
    
    prm.vax.strategy <- read.csv(file = file.vax.strategy, #  file.vax.strategy='prm_vax_strategy.csv' 
                                 as.is = T,
                                 strip.white = T) %>%
        create.list()
    
    
    # --- Merge all parameters:
    prm.hosp <- c(prm.hosp, prm.vax.features, prm.vax.strategy)
    
    # Only one single HCS:
    hcs <- list(prm.hosp)
    
    # Movements between HCS:
    mvt_hcs <- list(0) # no movements
    
    # --- Simulation parameters:
    sim.prm <- read.csv(file = file.prm.simul, 
                        as.is = TRUE) %>% 
        create.list()
    if(sim.prm$horizon_year > 0){
        sim.prm$horizon <- sim.prm$horizon_year * 365
        message(paste('Horizon conversion:',
                      sim.prm$horizon_year,
                      'years ->',
                      sim.prm$horizon,'days'))
    }
    
    # --- Disease parameters:
    disease.prm <- read.csv(file = file.disease, 
                            as.is = T, 
                            strip.white = T) %>% 
        create.list()
    
    # ---- Run Simulations ----
    
    message("Running simulations ...")
    n.mc    <- sim.prm$mc
    max.cpu <- parallel::detectCores()
    ncores  <- min(max.cpu, n.mc)
    sfInit(parallel = ncores>1, cpus = ncores)
    sfLibrary(abmcdiff, lib.loc = R.library.dir)
    sfExportAll()
    
    seedvec <- 1:n.mc
    
    res.list <- sfLapply(x = seedvec, 
                         fun = run_one, 
                         hcs = hcs, 
                         sim.prm = sim.prm, 
                         disease.prm = disease.prm, 
                         mvt_hcs = mvt_hcs)
    
    sfStop()
    
    tictac1 <- as.numeric(Sys.time())
    dt <- tictac1 - tictac
    message(paste(" ===> All simulations run in ", 
                  dt%/%60,'min ', 
                  round(dt%%60,0), "sec ==== "))
    
    
    # ---- Merge results -----
    
    message(paste('Merging',n.mc,'MC simulations ...'))
    res.mc <- lapply(X = 1:n.mc, FUN = digest_res, res.list)
    names(res.mc) <- paste0('mc',1:n.mc)
    df    <- do.call('rbind',res.mc)
    
    obj.to.save <- c('df')
    if(sim.prm$recordContacts){
        cd.mc  <- lapply(X = 1:n.mc, FUN = digest_res_cd, res.list)
        df.cd <- do.call('rbind',cd.mc)
        obj.to.save <- c(obj.to.save,'df.cd')
    }
    
    message('Merge done.')
    
    # Save for downstream use:
    message("Saving RData file ...")
    filename = paste0('simul-mc-run-simple.RData')
    save(list = obj.to.save, file = filename)
    message("RData file saved.")
    
    # ---- End stuff ----
    
    tictac2 <- as.numeric(Sys.time())
    dt <- tictac2 - tictac
    message(' |=====')
    message(paste(" |===== FULL EXECUTION RUN IN ", dt%/%60,'min ', 
                  round(dt%%60,0), "sec =====| "))
    message(' |=====')
    
    list.out <- list(df = df)
    if(sim.prm$recordContacts){
        list.out <- list(df = df, df.cd = df.cd)
    }
    
    return(list.out)
}

run_simulation_mc_new <- function(file.vax.strategy,
                              file.vax       = 'prm_vax.csv',
                              file.room.bed  = 'prm_room_bed.csv',
                              file.hosp.epi  = 'prm_hosp.csv',
                              file.prm.simul = 'prm_sim_comp.csv',
                              file.disease   = 'prm_disease.csv',
                              rngseed = 123456) {
    
    set.seed(rngseed)
    
    # ---- Set up ----
    # --- Health Care Settings definitions:
    file.room.hosp <- 'rooms_hosp.csv'
    
    z <- read.csv(file = file.room.bed ,#'prm_room_bed.csv', 
                  as.is = T,
                  strip.white = T)
    prm.room.bed <- create.list(z)
    
    r.hosp <- build_rooms(n.iso   = prm.room.bed$n.iso, 
                          n.wards = prm.room.bed$n.wards,
                          rooms.per.ward.mean = prm.room.bed$rooms.per.ward.mean, 
                          beds.per.room.mean  = prm.room.bed$beds.per.room.mean,
                          min.bed.per.room = prm.room.bed$min.bed.per.room,
                          max.bed.per.room = prm.room.bed$max.bed.per.room,
                          filename = file.room.hosp)
    summary_rooms(r.hosp)
    #plot_rooms(r.hosp)
    
    prm.hosp <- full.param.HCS(file_infrastruct = file.room.hosp, 
                               file_epi         = file.hosp.epi)
    
    # TO DO: CHANGE SINCE RANDOM GENERATION OF WARDS,ROOMS,BEDS.
    # REMOVE LIST ELEMENTS: room_type, ward, room_name, max_patient
    
    prm.hosp <- c(prm.hosp, prm.room.bed)
    
    # --- Vaccine features
    
    prm.vax <- create.list(read.csv(file = file.vax, #'prm_vax.csv', 
                                    as.is = T,
                                    strip.white = T))
    prm.vax.features <- list(vax_name = 'vaxificile',
                             vax_efficacy_lo = prm.vax$vax_efficacy_lo,
                             vax_efficacy_hi = prm.vax$vax_efficacy_hi,
                             vax_proba_prevent_colonization = prm.vax$vax_proba_prevent_colonization,
                             vax_halflife = prm.vax$vax_halflife)
    
    # --- Vaccination strategies
    
    prm.vax.strategy <- read.csv(file = file.vax.strategy, #  file.vax.strategy='prm_vax_strategy.csv' 
                                 as.is = T,
                                 strip.white = T) %>%
        create.list()
    
    
    # --- Merge all parameters:
    prm.hosp <- c(prm.hosp, prm.vax.features, prm.vax.strategy)
    
    # Only one single HCS:
    hcs <- list(prm.hosp)
    
    # Movements between HCS:
    mvt_hcs <- list(0) # no movements
    
    # --- Simulation parameters:
    sim.prm <- read.csv(file = file.prm.simul, 
                        as.is = TRUE) %>% 
        create.list()
    if(sim.prm$horizon_year > 0){
        sim.prm$horizon <- sim.prm$horizon_year * 365
    }
    
    # --- Disease parameters:
    disease.prm <- read.csv(file = file.disease, 
                            as.is = T, 
                            strip.white = T) %>% 
        create.list()
    
    # ---- Run Simulations ----
    
    message("Running simulations ...")
    n.mc    <- sim.prm$mc
    n.cpus  <- sim.prm$cpus
    n.cpus  <- ifelse(n.cpus==0, 
                      parallel::detectCores(), n.cpus)
    ncores  <- min(n.cpus, n.mc)
    
    sfInit(parallel = ncores>1, cpus = ncores)
    sfLibrary(abmcdiff, lib.loc = R.library.dir)
    sfExportAll()
    
    seedvec <- 1:n.mc
    
    res.list <- sfLapply(x = seedvec, 
                         fun = run_one, 
                         hcs = hcs, 
                         sim.prm = sim.prm, 
                         disease.prm = disease.prm, 
                         mvt_hcs = mvt_hcs)
    
    sfStop()
    
    tictac1 <- as.numeric(Sys.time())
    dt <- tictac1 - tictac
    message(paste(" ===> All simulations run in ", 
                  dt%/%60,'min ', 
                  round(dt%%60,0), "sec ==== "))
    
    
    # ---- Merge results -----
    
    message(paste('Merging',n.mc,'MC simulations ...'))
    res.mc <- lapply(X = 1:n.mc, FUN = digest_res, res.list)
    names(res.mc) <- paste0('mc',1:n.mc)
    df    <- do.call('rbind',res.mc)
    
    obj.to.save <- c('df')
    if(sim.prm$recordContacts){
        cd.mc  <- lapply(X = 1:n.mc, FUN = digest_res_cd, res.list)
        df.cd <- do.call('rbind',cd.mc)
        obj.to.save <- c(obj.to.save,'df.cd')
    }
    
    message('Merge done.')
    
    # Save for downstream use:
    message("Saving RData file ...")
    
    vax.tag <- strsplit(file.vax, split = '.', fixed = TRUE)[[1]][1]
    
    filename = paste0('simul-mc-',
                      prm.vax.strategy$vaxStrategy,
                      '-',
                      round(prm.vax.strategy$proba_vax_newAdmission*100),
                      '-',
                      vax.tag,
                      '.RData')
    save(list = obj.to.save, file = filename)
    message("RData file saved.")
    
    # ---- End stuff ----
    
    tictac2 <- as.numeric(Sys.time())
    dt <- tictac2 - tictac
    message(' |=====')
    message(paste(" |===== FULL EXECUTION RUN IN ", dt%/%60,'min ', 
                  round(dt%%60,0), "sec =====| "))
    message(' |=====')
    
    list.out <- list(df = df)
    if(sim.prm$recordContacts){
        list.out <- list(df = df, df.cd = df.cd)
    }
    
    return(list.out)
}



full.param.HCS <- function(file_infrastruct, file_epi) {
    prm0 <- read.csv(file_infrastruct, as.is = T)
    prm1 <- read.csv(file_epi, strip.white = T, as.is = T) %>% create.list()
    
    prm <- c(list(room_type   = prm0$type,
                  ward        = prm0$ward,
                  room_name   = prm0$name,
                  max_patient = prm0$max_patient),
             prm1)
    return(prm)
}

read_mvt_hcs <- function(fname) {
    x <- read.csv(fname, header=FALSE)
    # Return a list of rows:
    return(a <- as.list(as.data.frame(t(x))))
}

unit_df <- function(dt, hz, hn, mci) {
    
    thetime <- seq(dt,hz,by=dt)
    
    z <- data.frame(time = thetime ,
                    hcs_name = hn,
                    year = 1+thetime %/% 365,
                    mc = mci)
    return(z)
}

create_key <- function(df) {
    res <- df %>%
        mutate(key = paste(hcs_name,mc,year,time,sep='_'))
    return(res)
}

cumul_incidence <- function(df) {
    
    # retrieve simulation time step & horizon:
    dt <- unique(diff(df$time))
    dt <- dt[dt>0][1]
    hz <- max(df$time)
    
    df2 <- df %>%
        filter(type == 'patient') %>%
        filter(currentlySymptomatic) %>%
        mutate(year = 1 + time %/% 365) %>%
        group_by(time, hcs_name,  year, mc) %>%
        summarize(n = n_distinct(uid)) %>%
        arrange(hcs_name, mc, time) 
    
    # Names, year and monte-carlo:
    hnv <- unique(df2$hcs_name)
    mcv <- unique(df2$mc)

    # Create data frame with _all_ times
    # (not just the ones where incidence is positive)
    tmp <- list()
    m <- 1
    for(i in seq_along(hnv)){
        for(k in seq_along(mcv)){
            tmp[[m]] <- unit_df(dt,hz, hn=hnv[i], mci=mcv[k])
            m <- m+1
        }
    }
    df.full <- do.call('rbind', tmp) %>%
        create_key() %>%
        mutate(day.in.year = time %% 365) 
    
    df2  <- create_key(df2)
    df2a <- select(ungroup(df2), key,n)
    
    df.full2 <- left_join(df.full, df2a, by='key') 
    df.full2$n[is.na(df.full2$n)] <- 0
    
    df.cum <- df.full2 %>%
        arrange(hcs_name, year, mc, time) %>%
        group_by(hcs_name, year, mc) %>%
        mutate(cum.n = cumsum(n))
    
    return(df.cum)
}

#' Clean up the outputs (e.g., insert NAs, etc.)
cleanup_simul_df <- function(df) {
    
    dummy.number <- 999999
    
    # Insert NAs when patient has not been isolated:
    df$timeIsolationStart[df$timeIsolationStart <= -dummy.number] <- NA
    df$timeIsolationStart[df$timeIsolationEnd >= dummy.number] <- NA
    
    # Insert NAs when patient has not been symptomatic
    df$timeOnset[df$timeOnset <= -dummy.number] <- NA
    df$timeOnset[df$timeOnset >= dummy.number] <- NA
    
    return(df)
}


#' Digest simulation results.
digest_res <- function(i, res.list) {
    res <- res.list[[i]]
    res.pop   <- res[grepl('pop_',names(res))]
    
    # df.pop.dt <- lapply(X = res.pop, FUN = as.data.frame)  # <--- SLOWER!
    df.pop.dt <- lapply(X = res.pop, FUN = data.table::as.data.table) # <--- FASTER!
    
    # df.tmp <- do.call(what = 'rbind', args = df.pop.dt)  # <--- SLOWER!
    df.tmp <- data.frame(data.table::rbindlist(df.pop.dt)) # <--- FASTER!
    
    df.tmp$full_type <- paste(df.tmp$type, df.tmp$subtype, sep='_')
    df.tmp$mc <- i
    
    # Replace NA numbers (i.e. 9999999) by actual NAs:
    # Note: df.tmp[df.tmp==HUGENUMBER]  <- NA is _slower_.
    HUGENUMBER <- 999999
    df.tmp <- as.data.frame(lapply(df.tmp, function(x){
        replace(x, x == HUGENUMBER | x==-HUGENUMBER,NA)
    }))
    return(df.tmp)
}

#' Digest simulation results about contact durations.
digest_res_cd <- function(i, res.list){
    # Contact durations:
    res <- res.list[[i]]
    idx <- which(grepl('contact_duration',names(res)))
    cd.tmp <- list()
    for(j in 1:length(idx)){
        hcsname <- gsub('contact_duration_','',names(res[idx[j]]))
        cd.tmp[[j]] <-  data.frame(hcs_name = hcsname,
                                   cd = res[[ idx[j] ]], 
                                   mc = i)
    }
    # cd.mc[[i]] <- do.call('rbind',cd.tmp)
    z <- do.call('rbind',cd.tmp)
    return(z)
}

calc_patient_days_per_month <- function(df) {
    
    res <- df %>% 
        filter(type == 'patient') %>%
        mutate(month = 1 + time %/% 30) %>%
        group_by(mc, month, uid) %>%
        summarise(dos.mx = max(DoS)) %>%
        group_by(month, mc) %>%
        summarise(patient.days = sum(dos.mx))
    return(res)
}

calc_patient_days_per_year <- function(df) {
    
    res <- df %>% 
        filter(type == 'patient') %>%
        mutate(year = 1 + time %/% 365) %>%
        group_by(mc, year, uid) %>%
        summarise(dos.mx = max(DoS)) %>%
        group_by(year, mc) %>%
        summarise(patient.days = sum(dos.mx))
    return(res)
}

#' Check for non-sensical outputs
check_output <- function(df) {
    
    message('Checking outputs... ', appendLF = FALSE)
    
    df.hcw <- subset(df, type=='HCW')
    df.pat <- subset(df, type=='patient')

    # HCWs' attributes:
    stopifnot(sum(df.hcw$symptomatic)==0)
    stopifnot(sum(df.hcw$currentlySymptomatic)==0)
    #stopifnot(sum(df.hcw$totalIsolationDays, na.rm = TRUE)==0)
    stopifnot(sum(df.hcw$timeAdmission>0)==0)
    stopifnot(sum(!is.na(df.hcw$timeOnset))==0)
    stopifnot(sum(df.hcw$nRelapses>0)==0)
    
    # Times consistency:
    # idx <- !is.na(df$timeIsolationStart)
    # a <- df[idx,]
    # # stopifnot(sum(a$timeIsolationStart < a$timeOnset-1.1)==0)
    # 
    # idx <- !is.na(df$timeIsolationEnd)
    # a <- df[idx,]
    # stopifnot(sum(a$timeIsolationEnd < a$timeIsolationStart)==0)
    # 
    # 
    # # Relapses
    # idx <- which(df$nRelapses>0  & !df$symptomatic)
    # a <- df[idx,]
    # stopifnot(sum(is.na(df$timeRelapse[idx])) == 0)
    # stopifnot(sum(!df$symptomatic[idx]) == 0)
    
    # so far so good....
    message('OK.')
}


#' Calculate odds ratio and CI.
#' @param exp.ill Number exposed and ill
#' @param exp.well Number exposed and well
#' @param noexp.ill Number not exposed and ill
#' @param noexp.well Number not exposed and well
odds.ratio <- function(exp.ill,exp.well,noexp.ill,noexp.well, ci.level=0.95) {
    
    # exp.ill = 992
    # noexp.ill = 165
    # exp.well = 2260
    # noexp.well = 1017
    
    odds.ratio <- (exp.ill/noexp.ill) / (exp.well/noexp.well)
    
    se.log.or <- sqrt(1/exp.ill + 1/exp.well + 1/noexp.ill + 1/noexp.well)
    
    ci.lo.log <- log(odds.ratio) - 1.96*se.log.or
    ci.hi.log <- log(odds.ratio) + 1.96*se.log.or
    
    ci.hi <- exp(ci.hi.log)
    ci.lo <- exp(ci.lo.log)
    
    n.well <- exp.well+noexp.well
    n.ill <- exp.ill + noexp.ill
    return(c(or = odds.ratio, ci.lo=ci.lo, ci.hi=ci.hi,
             n.well = n.well, n.ill=n.ill))
}



clean_room_eff_distribution <- function(mean_clean_room_efficacy) {
    x <- seq(0,1,length.out = 1e3)
    m <- mean_clean_room_efficacy
    a <- 10
    b <- a*(1-m)/m
    y <- dbeta(x,a,b)
    plot(x,y, typ='l',lwd=5, main = 'Room cleaning: distribution of efficacy', 
         xlab='efficacy')
}

#' Monthly Incidence of symptomatic infections 
#' per 10,000 patient-days
incidence_patientdays_monthly <- function(df) {
    
    # Calculate patient-days for each month:
    df.pd <- calc_patient_days_per_month(df)%>%
        mutate(key = paste(mc,month,sep='_')) %>%
        ungroup()
    
    # Calculate number of cases per month: 
    df.monthly.inc <- df %>%
        filter(symptomatic & type == 'patient') %>%
        group_by(mc, uid) %>%
        summarise(ta = min(timeAcquisition, na.rm = TRUE)) %>%
        mutate(month = 1+ta%/%30) %>%
        arrange(mc, ta) %>%
        group_by(mc, month) %>%
        summarise(inc = n_distinct(uid))%>%
        mutate(key = paste(mc,month,sep='_')) %>%
        ungroup()
    
    dj <- left_join(df.pd %>%
                        select(key, patient.days),
                    df.monthly.inc, by='key') %>%
        mutate(inc.per.10000.pd = inc / patient.days * 10000)
    return(dj)
}


#' Annual Incidence of symptomatic infections 
#' per 10,000 patient-days
incidence_patientdays <- function(df) {
    
    # Calculate patient-days for each year:
    df.pd <- calc_patient_days_per_year(df)%>%
        mutate(key=paste(mc,year,sep='_')) %>%
        ungroup()
    
    # Calculate number of cases per year: 
    df.annual.inc <- df %>%
        filter(symptomatic & type == 'patient') %>%
        group_by(mc, uid) %>%
        summarise(ta = min(timeAcquisition, na.rm = TRUE)) %>%
        mutate(year = 1+ta%/%365) %>%
        arrange(mc, ta) %>%
        group_by(mc, year) %>%
        summarise(inc = n_distinct(uid))%>%
        mutate(key=paste(mc,year,sep='_')) %>%
        ungroup()
    
    dj <- left_join(df.pd %>%
                        select(key, patient.days),
                    df.annual.inc, by='key') %>%
        mutate(inc.per.10000.pd = inc / patient.days * 10000)
    return(dj)
}



prop_asymptom_dischrg <- function(z) {
    z <- z %>% 
        mutate(key = paste(vax.strategy, uid, mc, sep='_'))
    
    # Find the last time for each patient:
    df.tmax <- z %>%
        filter(type == 'patient') %>%
        group_by(vax.strategy, uid, mc ) %>%
        summarise(tmax = max(time)) %>%
        mutate(key = paste(vax.strategy, uid, mc, sep='_')) %>%
        ungroup() %>%
        select(key, tmax)
    
    # dataframe at discharge time:
    df.dischrg <- left_join(z,df.tmax, by='key') %>%
        filter(time == tmax)
    
    # calculate proportion asymptomatic at discharge time:
    df.pa <- df.dischrg %>%
        mutate(is.asympt = (isColonized & !symptomatic) ) %>%
        group_by(vax.strategy) %>%
        summarise(mn = mean(is.asympt),
                  sd = sd(is.asympt))
    
    return(df.pa)
    
}


ef1 <- function(x){
    # x %>% filter(!(!symptomatic & prdSymptom > 0))
    return(x[!(!x$symptomatic & x$prdSymptom > 0),])
}

ef2 <- function(x){x %>% filter(!(!symptomatic & totalIsolationDays > 0))}

ef3 <- function(x) {
    
    x2 = x[x$is.symptomatic, ]
    
    b = x2 %>%
        group_by(mc, uid) %>%
        summarize(rmx = max(nRelapses),
                  n = n()) 
    
    b$q   <- (b$n == b$rmx+1)
    b$key <- paste(b$mc,b$uid,sep='_')
    bf    <- b[b$q,]
    
    message(paste("\n ef3: proportion incomplete:",
                  round(mean(!b$q),6)))
    
    x$key <- paste(x$mc, x$uid,sep='_')
    
    return(x[x$key %in% bf$key,])
}

solidify <- function(x) { ef1(x) }

read_sensi <- function(file.names, vax.strategy.name) {
    # file.names = rdata.none
    # vax.strategy.name = 'none'
    n <- length(file.names)
    
    df.list <- list()
    
    for(i in 1:n){
        message(paste(i,'/',n,file.names[i] ))
        load(file.names[i])
        
        df.list[[i]] <- df %>%
            select(time, mc, uid, type, symptomatic, 
                   isColonized, acquisitionType) %>%
            mutate(vax.strategy = vax.strategy.name) %>%
            mutate(prm.val = prm.vec.k)
        rm(df)
    }
    res <- do.call('rbind', df.list)
    return(res)
}


response_sensi <- function(z) {
    
    z <- z %>% 
        mutate(key = paste(vax.strategy, uid, mc, prm.val, sep='_'))
    
    # Find the last time for each patient:
    df.tmax <- z %>%
        filter(type == 'patient') %>%
        group_by(vax.strategy, uid, mc, prm.val ) %>%
        summarise(tmax = max(time)) %>%
        mutate(key = paste(vax.strategy, uid, mc,prm.val, sep='_')) %>%
        ungroup() %>%
        select(key, tmax)
    
    # dataframe at discharge time:
    df.dischrg <- left_join(z,df.tmax, by='key') %>%
        filter(time == tmax)
    
    # calculate proportion asymptomatic at discharge time:
    df.pa <- df.dischrg %>%
        mutate(is.asympt = (isColonized & !symptomatic) ) %>%
        group_by(vax.strategy, prm.val) %>%
        summarise(mn = mean(is.asympt))
    
    # calculate proportion of environmental acquisition:
    df.env <- df.dischrg %>%
        mutate( envir.acq = (acquisitionType=='environment') ) %>%
        group_by(vax.strategy, prm.val) %>%
        summarise(mn = mean(envir.acq))
    df.env
    
    return(list(df.pa  = df.pa, 
                df.env = df.env))
}


