###
###   CREATE MAIN RESULTS FIGURES
###

source('utils.R')
library(dplyr);library(tidyr) ; library(ggplot2) #library(tidyverse) ; 
theme_set(theme_bw())

# ---- Helper Fcts ----

calc_inc_rdata <- function(rdatafile, qt.width = 80) {
    
    load(rdatafile)
    
    year.vax <- 5
    dj <- incidence_patientdays(df) 
    
    df.plot <- dj %>%
        mutate(year.after.vax = year - year.vax) %>%
        filter(year.after.vax >= 0) %>%
        group_by(year.after.vax) %>%
        summarise(m = mean(inc.per.10000.pd),
                  qlo = quantile(inc.per.10000.pd,probs = 0.5 - qt.width/200),
                  qhi = quantile(inc.per.10000.pd,probs = 0.5 + qt.width/200))
    
    df.plot$scen <- rdatafile
    return(df.plot)
}

calc_monthly_inc_rdata <- function(rdatafile, year.vax, qt.width = 80) {
    
    load(rdatafile)
    
    month.vax <- year.vax * 12
    dj <- incidence_patientdays_monthly(df) 
    
    df.plot <- dj %>%
        mutate(month.after.vax = month - month.vax) %>%
        filter(month.after.vax >= 0) %>%
        group_by(month.after.vax) %>%
        summarise(m = mean(inc.per.10000.pd),
                  qlo = quantile(inc.per.10000.pd,probs = 0.5 - qt.width/200),
                  qhi = quantile(inc.per.10000.pd,probs = 0.5 + qt.width/200))
    
    df.plot$scen <- rdatafile
    return(df.plot)
}

scenario_name <- function(x) {
    tmp <- gsub(pattern = 'simul-mc-',replacement = '', x = x)
    tmp <- gsub(pattern = '.RData',replacement = '', x = tmp)
    
    # To display in a particular order:
    idx0 <- grepl('none',tmp)
    idx1 <- grepl('frailty',tmp)
    tmp[idx0] <- paste0('0_',tmp[idx0])
    tmp[idx1] <- paste0('1_',tmp[idx1])
    
    return(tmp)
}

calc_sympt_prop <- function(datafile, 
                            year.min = 7) {
    message(paste('Loading',datafile,'...'),appendLF = F)
    load(datafile)
    message('done.')
    df$year <- df$time / 365
    
    a <- df %>%
        filter(year >= year.min ) %>%
        filter(nColonizations>0) %>%
        group_by(mc, acquisitionType) %>%
        summarise(p = mean(symptomatic)) %>%
        mutate(scen = scenario_name(datafile))
    return(a)
}

calc_sympt_prop2 <- function(datafile, 
                             year.min = 7) {
    message(paste('Loading',datafile,'...'),appendLF = F)
    load(datafile)
    message('done.')
    df$year <- df$time / 365
    
    a <- df %>%
        filter(year >= year.min ) %>%
        filter(type=='patient') %>%
        filter(symptomatic==TRUE) %>%
        group_by(mc,uid,acquisitionType) %>%
        summarise(n = n_distinct(uid)) %>%
        group_by(acquisitionType) %>%
        summarise(n.type = n()) %>%
        ungroup() %>%
        mutate(prop = n.type / sum(n.type)) %>%
        mutate(scen = scenario_name(datafile))
    return(a)
}

# ---- Figure Fcts ----


figure_month_inc_comp <- function(dfm) {
    g <- dfm %>%
        ggplot(aes(x = month.after.vax, 
                   y = m, 
                   colour = scenario, 
                   fill   = scenario))+
        geom_line(size=2)+
        geom_point(size=3, pch=21,fill='white', stroke=2)+ 
        geom_ribbon(aes(ymin=qlo, ymax=qhi), alpha=0.1, colour=NA)+
        ggtitle('Monthly incidence per 10,000 patient-days')+
        xlab('months after vaccination starts')+ 
        ylab('monthly inc / 10,000 patient-days')+
        scale_x_continuous(minor_breaks = FALSE)+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank())
    #g
    pdf('figure-monthly-inc-comp.pdf', width=10, height=6)
    plot(g)
    dev.off()
}

# s = "0_none-prm_vax_p_0"
# s = "1_frailty-prm_vax_p_5"

scen_readable <- function(s) {
    
    novax <- 'No vacc.'
    vax   <- 'Vacc.'
    n1    <- ifelse(grepl('none',s), novax , vax)
    n2    <- ''
    
    if(n1 == vax){
        v    <- 'frailty-'
        pos1 <- unlist(gregexpr(v,s))
        cvg  <- substr(s,pos1+nchar(v), pos1+nchar(v)+1)
        
        w   <- 'prm_vax_p_'
        pos <- unlist(gregexpr(w,s))
        val <- substr(s,pos+nchar(w), nchar(s))
        valpct <- paste0(as.numeric(val)*10, '%')
        
        n2 <- paste0('with cvg = ',cvg,'%; p = ',valpct)
    }
    res <- paste(n1, n2)
    return(res)
}


figure_1 <- function(df.plot,yr.rng) {
    
    df.plot$scen.plot <- NA
    for(i in 1:nrow(df.plot)){  #i=1
        df.plot$scen.plot[i] <- scen_readable(s = df.plot$scenario[i])
    }
    # Save data associated with figure:
    write.csv(df.plot, file = 'figure-1-data.csv',
              quote = F, row.names = F)
    
    g.inc.pd <- df.plot %>%
        ggplot(aes(x = year.after.vax, 
                   y = m, 
                   colour = scen.plot, 
                   fill = scen.plot))+
        geom_line(size=2)+
        geom_point(size=3, pch=21,fill='white', stroke=2)+ 
        geom_ribbon(aes(ymin=qlo, ymax=qhi), alpha=0.1, colour=NA)+
        ggtitle('FIGURE 1 - Annual incidence per 10,000 patient-days')+
        xlab('year after vaccination starts')+ 
        ylab('ann. inc / 10,000 patient-days')+
        scale_x_continuous(breaks=yr.rng[1]:yr.rng[2], minor_breaks = FALSE)+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank())
    
    pdf('figure-1.pdf', width=10, height=6)
    plot(g.inc.pd)
    dev.off()  
}

figure_2 <- function(df.plot) {
    
    df.plot2 <- df.plot %>%
        select(year.after.vax, m, scenario) %>%
        spread(scenario, m)
    
    idx0 <- which(grepl('none', names(df.plot2)))
    
    n <- ncol(df.plot2)
    rd.tag <- 'rel.diff.'
    for(i in 2:n){
        if(i != idx0){
            # print(i)
            tmp <- df.plot2[,i]/df.plot2[,idx0]-1
            df.plot2 <- cbind(df.plot2,tmp)
            names(df.plot2)[ncol(df.plot2)] <- paste0(rd.tag,names(df.plot2)[i])
        }
    }
    idx.rd <- which(grepl(rd.tag,names(df.plot2)))
    # Express in terms of 'reduction'
    df.plot2[,idx.rd] <- - df.plot2[,idx.rd] 
    # select only the columns needed:
    df.plot2 <- df.plot2[,c(1,idx.rd)]
    # Save data file associated with figure:
    write.csv(df.plot2, 
              file = 'figure-2-data.csv', 
              quote = F, row.names = F)
    
    df2 <- df.plot2 %>%
        gather('scen','reduction',2:ncol(df.plot2)) %>%
        mutate(scenario = gsub(rd.tag,'',scen))
    
    df2$scen.plot <- NA
    for(i in 1:nrow(df2)){
        df2$scen.plot[i] <- scen_readable(df2$scenario[i])
    }
    
    g2 <- ggplot(df2, 
                 aes(x = year.after.vax, 
                     y = reduction, 
                     colour = scen.plot))+
        geom_line(size=2)+
        geom_point(size=3, pch=21, fill='white', stroke=2)+
        xlab('Year after vaccination starts') + 
        ylab('Relative reduction')+
        ggtitle('FIGURE 2 - Relative reduction of symptomatic incidence vs. no vaccine')+
        guides(colour = guide_legend(title="Scenario"))+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank())
    pdf('figure-2.pdf', width=10, height=6)
    plot(g2)
    dev.off()
}

calc_rel_chng_propHAI <- function(a) {
    # Retrieve reference value:
    idx <- which(grepl('none', a$scen))
    ref.value <- a$p.type[idx]
    
    # Calculate relative difference:
    a$rel.diff <- a$p.type/ref.value - 1
    
    # Remove reference row:
    a <- filter(a, !grepl('none',scen))
    return(a)
}


figure_prop_HAI <- function(rdatafiles, 
                            ci.width, 
                            year.min) {
    
    tmp <- lapply(rdatafiles, calc_sympt_prop2, year.min)
    tmp2 <- do.call('rbind', tmp) %>%
        mutate(hai = ifelse(grepl('community',acquisitionType),
                            FALSE, TRUE))
    
    df <- tmp2 %>%
        group_by(scen, hai) %>%
        summarise(p.type = sum(prop)) %>%
        mutate(source = ifelse(hai,'Hospital', 'Community'))
    
    df$scen.plot <- NA
    for(i in 1:nrow(df)){
        df$scen.plot[i] <- scen_readable(df$scen[i])
    }
    
    
    g <- df %>%
        ggplot(aes(x=scen.plot, fill=source)) +
        geom_bar(aes(y=p.type), 
                 stat='identity',
                 position = 'dodge')+
        ggtitle('Source of symptomatic CDI')+
        xlab('Scenario')+ ylab('mean proportion') +
        theme(panel.grid.major.x = element_blank())+
        scale_fill_brewer(palette = 'Paired')
    
    pdf('figure-prop-hai.pdf', width=10, height=6)
    plot(g)
    dev.off()
    
}



figure_4 <- function(rdatafiles, ci.width, year.min) {
    
    tmp <- lapply(rdatafiles, calc_sympt_prop2, year.min)
    tmp2 <- do.call('rbind', tmp) %>%
        mutate(hai = ifelse(grepl('community',acquisitionType),
                            FALSE,TRUE))
    
    g4.a <- tmp2 %>%
        ggplot(aes(x=scen, fill=acquisitionType)) +
        geom_bar(aes(y=prop), 
                 stat='identity',
                 position = 'dodge')+
        ggtitle('Proportion of symptomatic HAI')+
        xlab('Scenario')+ ylab('mean proportion')
    
    df4b <- tmp2 %>%
        group_by(scen, hai) %>%
        summarise(p.type = sum(prop)) 
    
    g4.b <- df4b %>%
        ggplot(aes(x=scen, fill=hai)) +
        geom_bar(aes(y=p.type), 
                 stat='identity',
                 position = 'dodge')+
        ggtitle('Proportion of symptomatic HAI')+
        xlab('Scenario')+ ylab('mean proportion')
    
    p.hai <- calc_rel_chng_propHAI( filter(df4b, hai==TRUE))
    p.imp <- calc_rel_chng_propHAI( filter(df4b, hai==FALSE))
    df.reldiff.hai <- rbind(p.hai, p.imp)
    
    g4.c <- df.reldiff.hai %>%
        ggplot(aes(x = scen, y = rel.diff, colour= hai, group=hai))+
        geom_line(size=2)+ 
        geom_point(pch=21, size=3,stroke=2, fill='white')+
        geom_hline(yintercept = 0, linetype='dashed')+
        ylab('Relative difference')+
        ggtitle('Relative change of HAI and non-HAI proportions compared to No-Vaccine scenario')+
        scale_y_continuous(labels=scales::percent)
    g4.c
    
    pdf('figure-4a.pdf', width=10, height=6)
    plot(g4.a)
    plot(g4.b)
    plot(g4.c)
    dev.off()
    
}


# ---- RUN ----
do.run <- 1

if(do.run){
    # - Data
    rdatafiles <- system('ls simul-mc-*RData', intern = TRUE)
    ci.width   <- 95
    year.vax   <-  5
    
    df.list <- lapply(X = rdatafiles, 
                      FUN = calc_inc_rdata,
                      qt.width = ci.width)
    
    df.plot.all <- do.call('rbind',df.list)
    df.plot.all$scenario <- scenario_name(df.plot.all$scen)
    yr.rng <- range(df.plot.all$year.after.vax)
    
    # monthly incidence:
    df.list.m <- lapply(X = rdatafiles, 
                        FUN = calc_monthly_inc_rdata,
                        year.vax = year.vax,
                        qt.width = ci.width)
    dfm <- do.call('rbind',df.list.m)
    dfm$scenario <- scenario_name(dfm$scen)
    write.csv(dfm, 
              file='out-monthly-incidence-per-10kpatientdays.csv',
              quote = F, row.names = F)
    
    
    # - Figures
    
    figure_1(df.plot = df.plot.all, yr.rng)
    figure_2(df.plot = df.plot.all)
    try(figure_prop_HAI(rdatafiles, ci.width = ci.width, year.min = 7))
    
    figure_month_inc_comp(dfm)
    
    # figure_4(rdatafiles, ci.width = ci.width, year.min = 7)
    
    # Old stuff:
}





