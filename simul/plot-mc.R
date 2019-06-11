###
###  PLOT MONTE CARLO SIMULATIONS
###

library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
source('utils.R')

# ---- Preliminaries ----

t1 <- as.numeric(Sys.time())

# Load simulation output objects:
message('Loading data ...', appendLF = FALSE)
load('simul-mc-none-0-prm_vax_p_0.RData')
message(' done.')
df <- cleanup_simul_df(df)
check_output(df)

theme_set(theme_bw())

do.other.plots <- FALSE


# ---- Susceptibility index ----
g.si <- df %>%
    group_by(uid, type) %>%
    summarise(s = mean(sIdx)) %>%
    ggplot() +
    geom_histogram(aes(x=s, y=..density..,  
                       fill=type), 
                   alpha=0.8, bins=30)+
    ggtitle('Susceptibility Index') +
    xlab('Susceptibility index') +
    #coord_cartesian(xlim=c(0,1)) + 
    guides(fill=FALSE)+
    facet_grid(~type, scales = 'free')

g.ii <- df %>%
    group_by(uid, type) %>%
    summarise(s = mean(infectIdx)) %>%
    ggplot() +
    geom_histogram(aes(x=s, y=..density..,  
                       fill=type), 
                   alpha=0.8, bins=30)+
    ggtitle('Infection Index') +
    xlab('Infection index') +
    #coord_cartesian(xlim=c(0,1)) + 
    guides(fill=FALSE)+
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

df.incub <- df %>%
    filter(symptomatic==TRUE) %>%
    group_by(uid) %>%
    summarise(max.pi = max(prdIncubation)) 

m.incub  <- round(mean(df.incub$max.pi, na.rm=TRUE),1)
md.incub <- round(median(df.incub$max.pi, na.rm=TRUE),1)

g.incubPrd <- ggplot(df.incub) +
    geom_histogram(aes(x=max.pi),
                   fill='blue1', 
                   colour = 'blue3',
                   binwidth = 1) +
    ggtitle("Incubation period",
            paste('mean =',m.incub, '; median =',md.incub)) + 
    xlab('incubation period in days')
# g.incubPrd

# ---- Admission to onset ----
df.onset <- df %>%
    filter(!is.na(timeOnset)) %>%
    mutate(dt.adm.onset = timeOnset - timeAdmission)

# stats:
z <- df.onset %>%
    group_by(uid) %>%
    summarise(m = min(dt.adm.onset))

m.a2o  <- round(mean(z$m, na.rm = TRUE),1)
md.a2o <- round(median(z$m, na.rm = TRUE),1)

g.adm.onset <- NULL
if(nrow(df.onset)>0){
    g.adm.onset <- df.onset %>%
        group_by(uid)%>%
        summarize(max.dt = max(dt.adm.onset)) %>%
        ggplot()+
        geom_histogram(aes(x=max.dt), 
                       fill='lightgrey', colour='grey',
                       binwidth = 1)+
        ggtitle('Admission to onset interval',
                paste('mean =',m.a2o, '; median =',md.a2o))+
        xlab('admission to onset interval')+
        theme(panel.grid = element_blank())
}

# ---- Symptomatic period ----

df.sympt <- df %>%
    filter(symptomatic==TRUE) %>%
    group_by(uid) %>%
    summarise(max.ps = max(prdSymptom))

m.sympt  <- round(mean(df.sympt$max.ps, na.rm=TRUE),1)
md.sympt <- round(median(df.sympt$max.ps, na.rm=TRUE),1)
sd.sympt <- round(sd(df.sympt$max.ps, na.rm=TRUE),1)

g.symptPrd <-  ggplot(df.sympt) +
    geom_histogram(aes(x=max.ps),
                   fill='orange2', 
                   colour='orange3',
                   binwidth = 1) +
    ggtitle("Symptomatic period",
            paste('mean =',m.sympt, '; median =',md.sympt, '; sd =',sd.sympt)) +
    xlab('symptomatic period in days')+
    theme(panel.grid = element_blank())
# g.symptPrd

# ---- Contact durations ----

g.cd <- ggplot()+ggtitle('Contact durations NOT recorded')

if(exists('df.cd')){
    mcd <- df.cd %>%
        group_by(hcs_name) %>%
        summarise(m = mean(cd)*1440)
    means <- paste(c('means: ', 
                     paste(paste(mcd$hcs_name,'=',round(mcd$m,2)), 
                           collapse = ' ; ')),collapse='')
    
    g.cd <- ggplot(df.cd) + 
        geom_freqpoly(aes(x=cd * 1440, y=..density..), bins=50, size=2)+
        scale_x_log10(breaks=c(1,5,10,30,60,120))+ 
        scale_y_log10() +
        # geom_vline(xintercept = mcd, linetype=2) +
        ggtitle('HCW/patient contact durations',means) +
        xlab('Contact duration (minutes)')
}

# ---- Isolation period ----

g.isodurFirst <- df[df$symptomatic & df$isoDurFirst>0,] %>%
    ggplot()+
    geom_histogram(aes(x=isoDurFirst), binwidth = 1)+
    ggtitle('Isolation duration - 1st episode')

g.isodurRel <- df[df$symptomatic & df$isoDurRelapses>0,] %>%
    ggplot()+
    geom_histogram(aes(x=isoDurRelapses), binwidth = 1)+
    ggtitle('Isolation duration - Relapses')



# ---- Contamination index ----

g.contamIdx <- df %>%
    filter(isColonized) %>%
    group_by(mc, uid, symptomatic) %>%
    summarise(m = max(contamIdx)) %>%
    ggplot() +
    geom_density(aes(x=m, fill=symptomatic), colour=NA)+
    facet_wrap(~symptomatic, scales = 'free_y')+
    ggtitle('Maximum Contamination index among colonized')+
    xlab('Maximum contamIdx')


# ---- Relapses ----

g.rlp <- df %>%
    filter(symptomatic)%>%
    group_by(mc, uid) %>%
    summarise(nr = max(nRelapses)) %>%
    ggplot() + 
    geom_histogram(aes(x=nr, y=(..count..)/sum(..count..)), 
                   binwidth = 1 )+
    scale_y_continuous(breaks = seq(0,1,by=0.1))+
    xlab('Number of relapses')+
    ylab('frequency')+
    ggtitle('Number of relapses')

# ---- incidence -----

z = df %>%
    group_by(mc, uid) %>%
    summarise(bs = max(symptomatic))
message(paste('Full df> cases per 10,000 patients: ', mean(z$bs) * 10000))

# symptomatic incidence per patient days
dj <- incidence_patientdays(df)
g.inc.pd <- dj %>%
    ggplot(aes(x=factor(year), y=inc.per.10000.pd))+
    geom_boxplot(fill='lightgrey')+
    ggtitle('Annual incidence per 10,000 patient-days')+
    xlab('year')+ ylab('ann. inc / 10,000 patient-days')


if(do.other.plots){
    
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
    
    message("Plot object population done.")
    
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
    
    message("Plot object traffic done.")
    
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
    
    
    
    
    
    # by month to double checK:
    djm <- incidence_patientdays_monthly(df)
    g.inc.pd.m <- djm %>%
        ggplot(aes(x=factor(month), y=inc.per.10000.pd))+
        geom_boxplot(fill='lightgrey')+
        ggtitle('Monthly incidence per 10,000 patient-days')+
        xlab('month')+ ylab('monthly inc / 10,000 patient-days')
    
    
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
    
    
    # ---- Prop Symptomatic ----
    
    
    tmp <-  df %>%
        filter(type == 'patient') %>%
        group_by(mc, uid) %>%
        summarise(s = max(symptomatic))%>%
        group_by(mc) %>%
        summarise(ms = mean(s))
    
    
    tmp <-  df %>%
        filter(type == 'patient') %>%
        group_by(mc, uid, symptomatic) %>%
        summarise(col = max(isColonized))%>%
        filter(col==1) %>%
        group_by(mc, uid) %>%
        summarise(s = max(symptomatic))%>%
        group_by(mc) %>%
        summarise(ms = mean(s))
    
    
    g.ps <- ggplot(tmp)+
        geom_boxplot(aes(x=factor(1), y=ms), fill='lightgrey')+
        ggtitle('Symptomatic Proportion')+
        theme(axis.text.x  = element_blank(),
              axis.title = element_blank())
    
    g.pts <- df %>%
        filter(type == 'patient') %>%
        group_by(mc, uid) %>%
        summarise(rho = max(prdSymptom)/max(DoS)) %>%
        filter(rho>0) %>%
        ggplot() + 
        geom_boxplot(aes(x=factor(1), y=rho))
    
    # pts <- mean(a$rho[a$rho>0])
    
    # ---- Acquisition type ----
    
    z <- df %>%
        group_by(mc, uid, acquisitionType) %>%
        summarise(n = n_distinct(acquisitionType)) %>%
        filter(n==1) %>%
        group_by(acquisitionType, mc) %>%
        summarise(nat = sum(n))
    
    g.acqt <- z %>% 
        filter(acquisitionType != "NA") %>%
        group_by(acquisitionType) %>%
        summarise(y = mean(nat),
                  yy = sd(nat)) %>%
        ggplot()+
        geom_pointrange(aes(x=acquisitionType, 
                            y=y, 
                            ymin = y - yy,
                            ymax = y + yy,
                            colour=acquisitionType),
                        size=3)+
        ggtitle('Mean (+/- sd) number of acquisition types')+
        guides(colour=FALSE)+
        theme(text = element_text(size = 20))+
        xlab('')+ylab('')
    
    
    df.acqtyp <- df %>%
        filter(type=='patient') %>%
        filter(symptomatic==TRUE) %>%
        group_by(mc,uid,acquisitionType) %>%
        summarise(n = n_distinct(uid)) %>%
        group_by(acquisitionType) %>%
        summarise(n.type = n()) %>%
        ungroup() %>%
        mutate(prop = n.type / sum(n.type))
    
    g.acqt.sympt <- df.acqtyp %>%
        ggplot(aes(x=acquisitionType, y=n.type))+
        geom_bar(stat='identity', fill='grey80')+
        geom_text(aes(label=n.type), size=8)+
        ggtitle('Symptomatic infections by source', 'raw numbers')
    
    g.acqt.sympt.prop <- df.acqtyp %>%
        ggplot(aes(x=acquisitionType, y=prop))+
        geom_bar(stat='identity', fill='grey80')+
        geom_text(aes(label=round(prop,2)), size=8)+
        ggtitle('Symptomatic infections by source', 'proportion')
    
    message('Plot object acquisition type done.')
    
    # ---- Contamination ----
    
    # Rooms:
    
    tb <- 7 # time bucket (else plot has many points ==> slow to display)
    
    rc <- df %>%
        filter(current_room_uid < 9999) %>%
        filter(mc == 1 ) %>%
        mutate(timeb = round(time/tb)*tb ) %>%
        group_by(timeb, current_room_uid) %>%
        summarise(y = mean(current_room_contamIdx))
    
    g.cont <- ggplot(rc, aes(x=timeb, y=y)) +
        geom_line(size=0.5, colour='gray') + 
        geom_point(size=0.5) +
        facet_wrap(~current_room_uid) +
        ggtitle('Room environmental contamination (MC=1)') + 
        xlab('bucketed time')+
        ylab('Load')
    
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
    
    g.oc <- ggplot(oc, aes(x=time, y=y)) +
        geom_line(size=2) +
        geom_line(data = mpr, aes(x=time, y=el), colour='blue', linetype=3)+
        ggtitle('Other environmental contamination')
    
    # ---- Shedding index ----
    
    if(FALSE){
        ts.shi <- df %>%
            group_by(time, uid, type) %>%
            summarise(shi = mean(contamIdx))
        
        
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
        
        g.vax.dos <- df %>%
            filter(isColonized==TRUE) %>%
            filter(type=='patient') %>%
            ggplot(aes(y=DoS, x=isVax, fill=isVax))+
            geom_boxplot(size=1)+
            ggtitle('Vaccination impact on DoS','Among all colonized patients')
    }
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

# ---- PLOT ALL ----

message("Saving plot to `plot-mc.pdf` ...", appendLF = FALSE)

pdf('plot-mc.pdf', width = 15, height = 10)

try({
    grid.arrange(g.si, g.ii, g.dos,
                 g.symptPrd, g.incubPrd, 
                 g.adm.onset, g.contamIdx, 
                 g.cd, g.rlp, g.isodurFirst, g.isodurRel)
})

plot_cum_prev_year(df, colnz.or.sympt = 'symptomatic') 
plot(g.inc.pd)


if(do.other.plots){
    try(plot(g.pop))
    
    try(grid.arrange(g.patientdays, g.adm, nrow=1))
    
    # try({
    #     grid.arrange(plot_cum_slope(df, colnz.or.sympt = 'colonized'),
    #              plot_cum_slope(df, colnz.or.sympt = 'symptomatic'),
    #              nrow=1)
    # })
    
    plot_cum_prev_year(df, colnz.or.sympt = 'colonized') 
    plot_cum_prev_year(df, colnz.or.sympt = 'symptomatic') 
    
    plot(g.inc.pd.m)
    plot(g.cr)
    
    plot(g.ps)
    plot(g.pts)
    
    plot(g.acqt)
    grid.arrange(g.acqt.sympt, 
                 g.acqt.sympt.prop, 
                 nrow=1)
    
    plot(g.cont)
    plot(g.oc)
    
    if(length(uv)==2){
        plot(g.vax.t)
        grid.arrange(try(plot_oddsratio(or.col, 'Colonization')),
                     try(plot_oddsratio(or.sympt, 'Symptomatic')),
                     ncol=1)
        plot(g.vax.dos)
    }
}
message('done.')

dev.off()


# ---- End ----
t2 <- as.numeric(Sys.time())
dt <- round( (t2-t1)/60, 1)
msg <- paste('--- ALL PLOTS DONE in',dt, 'minutes ---')
message('---')
message(msg)

