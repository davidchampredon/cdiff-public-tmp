###
###  PLOT MONTE CARLO SIMULATIONS
###

library(tidyverse)
library(gridExtra)
source('utils.R')


plot_one_strategy <- function(rdatafile) {
    
    # ---- Preliminaries ----
    
    t1 <- as.numeric(Sys.time())
    
    # Load simulation output objects:
    message(paste('Loading data',rdatafile,'...'), appendLF = FALSE)
    load(file = rdatafile)
    message(' done.')
    
    df <- cleanup_simul_df(df)
    check_output(df)
    
    # Does not deal with more than 1 HCS:
    stopifnot(length(unique(df$hcs_name))==1)
    
    theme_set(theme_bw())
    
    vax.strategy.name <- gsub(pattern = 'simul-mc-', replacement = '', 
                              x = rdatafile, fixed = TRUE)
    vax.strategy.name <- gsub(pattern = '.RData', replacement = '', 
                              x = vax.strategy.name, fixed = TRUE)
    
    file.plot <- paste0('plot-mc-',vax.strategy.name,'.pdf')
    pdf(file.plot, width = 15, height = 10)
    
    
    # ---- Susceptibility index ----
    g.si <- df %>%
        group_by(uid, type, hcs_name) %>%
        summarise(s = mean(sIdx)) %>%
        ggplot() +
        geom_histogram(aes(x=s, y=..density..,  fill=type), 
                       alpha=0.8, binwidth = 0.05)+
        ggtitle('Susceptibility Index') +
        xlab('Susceptibility index') +
        coord_cartesian(xlim=c(0,1)) + guides(fill=FALSE)+
        facet_grid(~type, scales = 'free')
    
    # ---- DoS ----
    
    df.dos <- df %>%
        filter(type=='patient') %>%
        select(DoS,symptomatic)
    
    mean.dos <- df.dos %>%
        group_by(symptomatic) %>%
        summarise(m = mean(DoS),
                  md = median(DoS))
    
    g.dos <- df.dos %>%
        ggplot() +
        geom_density(aes(x=DoS,fill=symptomatic), 
                     colour=FALSE, 
                     alpha=0.5)+
        geom_vline(data=mean.dos, 
                   aes(xintercept=m, colour=symptomatic),
                   size=2, alpha=0.5)+
        geom_point(data=mean.dos, 
                   aes(x=md,y=0, colour=symptomatic),
                   size=2, alpha=0.5)+
        ggtitle(paste("Duration of stay - mean =",
                      paste(round(mean.dos$m,1),collapse = ' & '))) + 
        xlab('Duration of stay in days')+
        scale_x_log10(breaks=c(1,2,5,7,10,14,30,60,90), minor_breaks=FALSE)+
        theme(panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank())
    # g.dos
    
    
    # ---- Incubation ----
    g.incubPrd <- df %>%
        filter(symptomatic==TRUE) %>%
        group_by(uid) %>%
        summarise(max.pi = max(prdIncubation)) %>%
        ggplot() +
        geom_histogram(aes(x=max.pi),
                       fill='blue1', 
                       colour = 'blue3',
                       binwidth = 1) +
        ggtitle("Incubation period") + 
        xlab('incubation period in days')
    
    
    # ---- Admission to onset ----
    df.onset <- df[!is.na(df$timeOnset),]
    
    g.adm.onset <- NULL
    if(nrow(df.onset)>0){
        g.adm.onset <- df.onset %>%
            mutate(dt.adm.onset = timeOnset - timeAdmission) %>%
            group_by(uid)%>%
            summarize(max.dt = max(dt.adm.onset)) %>%
            ggplot()+
            geom_histogram(aes(x=max.dt), 
                           fill='lightgrey', colour='grey',
                           binwidth = 1)+
            ggtitle('Admission to onset interval')+
            theme(panel.grid = element_blank())
    }
    
    # ---- Symptomatic period ----
    
    g.symptPrd <- df %>%
        filter(symptomatic==TRUE) %>%
        group_by(uid) %>%
        summarise(max.ps = max(prdSymptom)) %>%
        ggplot() +
        geom_histogram(aes(x=max.ps),
                       fill='orange2', 
                       colour='orange3',
                       binwidth = 1) +
        ggtitle("Symptomatic period") +
        xlab('symptomatic period in days')+
        theme(panel.grid = element_blank())
    
    # ---- Contact durations ----
    
    mcd <- df.cd %>%
        group_by(hcs_name) %>%
        summarise(m = mean(cd)*1440)
    means <- paste(c('means: ', 
                     paste(paste(mcd$hcs_name,'=',round(mcd$m,2)), collapse = ' ; ')),collapse='')
    
    g.cd <- ggplot(df.cd) + 
        geom_freqpoly(aes(x=cd * 1440, y=..density..), bins=50, size=2)+
        scale_x_log10(breaks=c(1,5,10,30,60,120))+ 
        scale_y_log10() +
        # geom_vline(xintercept = mcd, linetype=2) +
        ggtitle('HCW/patient contact durations',means) +
        xlab('Contact duration (minutes)')
    
    
    # ---- Isolation period ----
    
    g.tot.iso <- df %>%
        filter(totalIsolationDays>0) %>%
        group_by(uid) %>%
        summarise(x = max(totalIsolationDays)) %>%
        ggplot()+
        geom_histogram(aes(x=x), binwidth = 1)+
        ggtitle('Total isolation duration')+
        xlab('isolation days')
    
    g.onset.iso <-  df %>%
        filter(timeIsolationStart>0) %>%
        mutate(dt = round(timeIsolationStart - timeOnset)) %>%
        group_by(uid) %>%
        summarize(max.dt = max(dt)) %>%
        ggplot()+
        geom_histogram(aes(x=max.dt), binwidth = 1) +
        scale_x_continuous(minor_breaks = FALSE)+
        ggtitle('Onset to isolation lag')
    
    
    # ---- Plots in grid ----
    grid.arrange(g.si, g.dos,
                 g.symptPrd, g.incubPrd, 
                 g.adm.onset, g.onset.iso, 
                 g.tot.iso, g.cd)
    
    
    # ---- Population ----
    
    g.pop <- df %>%
        filter(time>10) %>%
        group_by(mc, time, type) %>%
        summarise(n = n_distinct(uid)) %>%
        group_by(time,type) %>%
        summarise(m  = mean(n)) %>%
        ggplot(aes(x=time)) + 
        geom_line(aes(y=m), size=1)+
        facet_wrap(~type, scale='free')+
        ggtitle('Population time series (mean)')
    plot(g.pop)
    
    # ---- Traffic ----
    
    df.adm <- df %>%
        group_by(timeAdmission,mc)%>%
        summarise(n = n_distinct(uid)) %>%
        mutate(yr = 1+timeAdmission%/%365) %>%
        group_by(yr,mc) %>%
        summarise(adm.yr = sum(n))
    
    g.adm <- df.adm %>%
        ggplot(aes(x=factor(yr), y=adm.yr)) +
        geom_boxplot(fill='grey')+
        ggtitle('Annual number of admissions')+
        xlab('year') + ylab('number of admissions')
    
    
    g.patientdays <- df %>%
        calc_patient_days_per_year() %>%
        ggplot(aes(x=factor(year), y=patient.days)) +
        geom_boxplot(fill='lightgrey')+
        xlab('year') +
        ggtitle('Patient days')
    
    grid.arrange(g.patientdays, g.adm, nrow=1)
    
    
    # ---- Time series of prevalence ----
    
    linreg <- function(x,y) {
        z <- lm(formula = y ~ x + 0)
        k <- coef(summary(z))[1]
        return(k)
    }
    
    plot_cum_slope <- function(df, colnz.or.sympt) {
        # colnz.or.sympt='symptomatic'
        if(grepl('col',colnz.or.sympt)) df2 <- filter(df,isColonized==TRUE)
        if(grepl('sym',colnz.or.sympt)) df2 <- filter(df,symptomatic==TRUE)
        
        df3 <- df2 %>%
            group_by(mc, uid, full_type) %>%
            summarise(ta = min(timeAcquisition, na.rm = TRUE)) %>%
            mutate(year = 1+ta%/%365) %>%
            mutate(day.of.year = ta%%365) %>%
            arrange(mc, ta) %>%
            mutate(dummy = 1) %>%
            group_by(mc,year, full_type)%>%
            mutate(n.cum = cumsum(dummy)) 
        
        df4 <- df3 %>%
            group_by(mc, year) %>%
            summarise(slope.cum.inc = linreg(day.of.year, n.cum))
        
        g1 <- df4  %>%
            ggplot()+
            geom_density(aes(x=slope.cum.inc), fill='tomato1',alpha=.5)+
            geom_vline(xintercept = mean(df4$slope.cum.inc), linetype='dashed',size=2)+
            ggtitle(paste('Slope of cumulative incidence - ',
                          ifelse(grepl('col',colnz.or.sympt),'colonized','symptomatic')),
                    'Across all years and MCs')
        
        return(g1)
    }
    
    plot_cum_prev_year <- function(df, colnz.or.sympt) {
        # colnz.or.sympt='symptomatic'
        if(grepl('col',colnz.or.sympt)) df2 <- filter(df,isColonized==TRUE)
        if(grepl('sym',colnz.or.sympt)) df2 <- filter(df,symptomatic==TRUE)
        
        df3 <- df2 %>%
            group_by(mc, uid, full_type) %>%
            summarise(ta = min(timeAcquisition, na.rm = TRUE)) %>%
            mutate(year = 1+ta%/%365) %>%
            mutate(day.of.year = ta%%365) %>%
            arrange(mc, ta) %>%
            mutate(dummy = 1) %>%
            group_by(mc,year, full_type)%>%
            mutate(n.cum = cumsum(dummy)) 
        
        g <- df3 %>%
            ggplot(aes(x=day.of.year, y=n.cum, colour=factor(year)))+
            geom_step(size=2,alpha=.8)+
            facet_grid(full_type~mc)+
            scale_y_continuous(minor_breaks = FALSE)+
            ggtitle(paste('Cumulative incidence -',colnz.or.sympt))+
            xlab('Acquisition time')+ylab('Cumulative incidence') +
            scale_color_brewer(palette = 'Blues')
        plot(g)
    }
    
    plot_cum_prev_vax <- function(df, colnz.or.sympt) {
        # colnz.or.sympt='symptomatic'
        if(grepl('col',colnz.or.sympt)) df2 <- filter(df,isColonized==TRUE)
        if(grepl('sym',colnz.or.sympt)) df2 <- filter(df,symptomatic==TRUE)
        
        df3 <- df2 %>%
            group_by(mc, uid, full_type,isVax) %>%
            summarise(ta = min(timeAcquisition, na.rm = TRUE)) %>%
            arrange(mc,isVax, ta) %>%
            mutate(dummy = 1) %>%
            group_by(mc, isVax)%>%
            mutate(n.cum = cumsum(dummy))%>%
            mutate(vax.status = ifelse(isVax,'Vaccinated', 'No vaccine'))
        
        g <- df3 %>%
            ggplot(aes(x=ta, y=n.cum, size=full_type, colour=vax.status))+
            geom_step(size=2,alpha=.8)+
            facet_wrap(~mc)+
            scale_y_continuous(minor_breaks = FALSE)+
            ggtitle(paste('Cumulative incidence -',colnz.or.sympt))+
            xlab('Acquisition time')+ylab('Cumulative incidence')
        plot(g)
    }
    
    grid.arrange(plot_cum_slope(df, colnz.or.sympt = 'colonized'),
                 plot_cum_slope(df, colnz.or.sympt = 'symptomatic'),
                 nrow=1)
    
    plot_cum_prev_vax(df, colnz.or.sympt = 'colonized') 
    plot_cum_prev_vax(df, colnz.or.sympt = 'symptomatic') 
    plot_cum_prev_year(df, colnz.or.sympt = 'colonized') 
    plot_cum_prev_year(df, colnz.or.sympt = 'symptomatic') 
    
    # by room
    ts.c.r <- df %>%
        filter(mc==1) %>%
        filter(type == 'patient') %>%
        group_by(time, current_room_uid) %>%
        summarise(n = length(uid),
                  nc = sum(isColonized),
                  ns = sum(symptomatic))
    g.cr <- ggplot(ts.c.r, aes(x=time)) +
        geom_step(aes(y=n),size=0.5, colour='black') +
        geom_step(aes(y=nc),size=1.2, colour='orange') +
        geom_step(aes(y=ns),size=1.2, colour='red') +
        facet_wrap(~current_room_uid, scales = 'fixed') +
        ggtitle('Patient occupancy, Colonized & Symptomatic per room (mc=1)')
    plot(g.cr)
    
    
    # ---- Contamination ----
    
    # Rooms:
    
    tb <- 1 # time bucket (else plot has many points ==> slow to display)
    
    # rc <- df %>%
    #     filter(current_room_uid < 9999) %>%
    #     filter(mc == 1 ) %>%
    #     mutate(timebckt = round(time/tb)*tb) %>%
    #     group_by(timebckt, current_room_uid) %>%
    #     summarise(el = mean(current_room_contamIdx),
    #               n = length(current_room_contamIdx))
    
    # ggplot(rc, aes(x=time, y=el)) +
    #     geom_line(size=0.5, colour='gray') + geom_point(size=0.5) +
    #     facet_wrap(~current_room_uid) +
    #     ggtitle('Room environmental contamination') + 
    #     xlab('time')+
    #     ylab('Load')
    
    rc <- df %>%
        filter(current_room_uid < 9999) %>%
        filter(mc == 1 ) 
    
    ggplot(rc, aes(x=time, y=current_room_contamIdx)) +
        geom_line(size=0.5, colour='gray') + 
        geom_point(size=0.5) +
        facet_wrap(~current_room_uid) +
        ggtitle('Room environmental contamination') + 
        xlab('time')+
        ylab('Load')
    
    unique(df$current_room_uid)
    
    # Individuals:
    if(length(unique(df$uid))<100){
        ici <- df %>%
            group_by(full_type, uid, time) %>%
            summarize(idx = mean(contamIdx),
                      col = mean(isColonized),
                      sym = mean(currentlySymptomatic),
                      n = length(contamIdx))
        
        ici$status <- 'carrier'
        ici$status[ici$col==1] <- 'colonized'
        ici$status[ici$sym==1] <- 'symptomatic'
        
        ggplot(ici, aes(x=time, y=idx,
                        colour=factor(status))) +
            #geom_line(size=0.5) +
            geom_point()+
            facet_wrap(full_type~uid)
    }
    
    # Environmental:
    mpr <- df %>%
        filter(type=='patient') %>%
        group_by(time) %>%
        summarize(el = mean(current_room_contamIdx),
                  n = length(current_room_contamIdx))
    
    oc <- df %>%
        group_by(time) %>%
        summarize(y = mean(contamIdx_other))
    ggplot(oc, aes(x=time, y=y)) +
        geom_line(size=2) +
        geom_line(data = mpr, aes(x=time, y=el), colour='blue', linetype=3)+
        ggtitle('Other environmental contamination')
    
    # ---- Shedding index ----
    ts.shi <- df %>%
        group_by(time, uid, type) %>%
        summarise(shi = mean(contamIdx))
    
    if(FALSE){
        ggplot(ts.shi,aes(x=time, y=shi, colour=factor(type))) +
            geom_line() + #geom_point()+
            ggtitle('contamIdx') +
            facet_wrap(~uid, scales='free_y') +
            guides(colour=FALSE)
    }
    
    # ---- Vaccination ----
    
    # Check is there is at least one vaccinated individual:
    uv <- unique(df$isVax)
    
    if(length(uv)==2){
        df.vax.t <- df %>%
            filter(type == 'patient') %>%
            group_by(mc, time, isVax) %>%
            summarise(n = n_distinct(uid)) %>%
            group_by(time,isVax) %>%
            summarize(m = mean(n)) %>%
            spread(isVax,m) %>%
            mutate(p.vax = `TRUE`/(`TRUE`+`FALSE`))
        
        g.vax.t <- df.vax.t %>%
            ggplot() +
            geom_line(aes(x=time, y=p.vax))+
            coord_cartesian(ylim=c(0,1))+
            scale_y_continuous(breaks = seq(0,1,by=0.2))+
            geom_hline(aes(yintercept = mean(p.vax)),linetype=2)+
            ggtitle('Mean vaccination proportion (across all MCs)')
        plot(g.vax.t)
        
        vax.odds.col <- df %>%
            group_by(mc,isVax, isColonized) %>%
            summarise(n = n_distinct(uid)) 
        
        vax.odds.sympt <- df %>%
            group_by(mc,isVax, symptomatic) %>%
            summarise(n = n_distinct(uid)) 
        
        calc_or <- function(i, df.odds, varname) {
            dfi <- filter(df.odds, mc==i)
            
            exp.ill    <- dfi$n[dfi$isVax==TRUE  & dfi[,varname]==TRUE]
            exp.well   <- dfi$n[dfi$isVax==TRUE  & dfi[,varname]==FALSE]
            noexp.ill  <- dfi$n[dfi$isVax==FALSE & dfi[,varname]==TRUE]
            noexp.well <- dfi$n[dfi$isVax==FALSE & dfi[,varname]==FALSE]
            
            or <- odds.ratio(exp.ill, exp.well,noexp.ill,noexp.well)
            nom <- names(or)
            check <- sum(sapply(c(exp.ill, exp.well,noexp.ill,noexp.well),length))
            if(check<4) or <- c(or=NA,ci.lo=NA,ci.hi=NA,n.well=NA,n.ill=NA)
            return(or)
        }
        
        or.col   <- sapply(X=unique(df$mc), FUN = calc_or, df.odds=vax.odds.col, varname='isColonized') %>% 
            t() %>%
            as.data.frame()
        or.sympt <- sapply(X=unique(df$mc), FUN = calc_or, df.odds=vax.odds.sympt, varname='symptomatic') %>% 
            t() %>%
            as.data.frame()
        
        plot_oddsratio <- function(df.oddsratio,title='') {
            g <- df.oddsratio %>%
                mutate(mc=1:nrow(df.oddsratio)) %>%
                ggplot(aes(x=factor(mc)))+
                geom_hline(yintercept = 1, 
                           linetype='dashed', colour='red',size=2)+
                geom_segment(aes(xend=factor(mc),y=ci.lo,yend=ci.hi),
                             size=1.3)+
                geom_point(aes(y=or, size=n.ill/n.well),
                           pch=21,fill='lightgrey',stroke=1)+
                guides(size=FALSE)+
                ggtitle(paste('Odds ratio when vaccinated -',title))+
                ylab('Odds ratio (CI)')+
                scale_y_log10()+
                theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank()
                )
            return(g)
        }
        grid.arrange(try(plot_oddsratio(or.col, 'Colonization')),
                     try(plot_oddsratio(or.sympt, 'Symptomatic')),
                     ncol=1)
        
        g.vax.dos <- df %>%
            filter(isColonized==TRUE) %>%
            filter(type=='patient') %>%
            ggplot(aes(y=DoS, x=isVax, fill=isVax))+
            geom_boxplot(size=1)+
            ggtitle('Vaccination impact on DoS','Among all colonized patients')
        plot(g.vax.dos)
    }
    
    
    
    # ---- End ----
    
    dev.off()
    t2 <- as.numeric(Sys.time())
    dt <- round( (t2-t1)/60, 1)
    msg <- paste('--- Plots for',rdatafile,'done in',dt, 'minutes ---')
    message('---')
    message(msg)
    
    
    
    
}

plot_all_strategy <- function() {
    rdata.list <- system('ls simul-mc-*.RData', intern = TRUE)
    x <- lapply(X = rdata.list, FUN = plot_one_strategy)
}


run.test <- TRUE

if(run.test){
    plot_all_strategy()
}

