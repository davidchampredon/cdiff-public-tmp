library(dplyr)
library(tidyr) 
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)

source('utils.R')

# ---- Functions ----
total_memory_used <- function() {
    
}


load_all_simul <- function(svax) {
    
    # Retrieve all simulations:
    rdata.list <- system('ls simul-mc-*.RData', intern = TRUE)
    
    # Extract the ones selected:
    idx <- sapply(svax, grep, rdata.list)
    rdata.list <- rdata.list[idx]
    print(rdata.list)
    
    vaxstrgy <- gsub(pattern = 'simul-mc-', replacement = '', 
                     x = rdata.list, fixed = TRUE)
    vaxstrgy <- gsub(pattern = '.RData',    replacement = '', 
                     x = vaxstrgy, fixed = TRUE)
    
    df.list <- list()
    n <- length(rdata.list)
    if(n<2){
        msg <- 'Vaccination strategies selected found in less than 2 RData files: Aborting!'
        message(msg); print(msg); stop()
    }
    for(i in 1:n){ #i=1
        print(paste(i,'/',n,': Loading',rdata.list[i]))
        load(rdata.list[i])
        df$vax.strategy <- vaxstrgy[i]
        df.list[[i]] <- df
        rm(df)
        mem <- format(object.size(df.list),units = 'MB')
        print(paste('Memory df.list:',mem))
    }
    # Too slow: df.all <- do.call('rbind',df.list) 
    # Use data.table's "rbindlist" instead
   df.all <- data.table::rbindlist(df.list)
   
   return(df.all)
}

plot_tmp <- function(z) {
    g <- z %>%
        filter(mc <= 9) %>%
        group_by(vax.strategy, mc, time) %>%
        summarize(n = sum(symptomatic)) %>%
        ggplot(aes(x=time, y=n, colour = vax.strategy)) + 
        geom_line(size=1, alpha=0.6)+
        facet_grid(vax.strategy~mc)
    plot(g)
    
    tmp <- z %>%
        group_by(vax.strategy, mc, time) %>%
        summarize(n = sum(symptomatic)) %>%
        spread(key = vax.strategy, value = n)
    nc <- ncol(tmp)
    
    tmp <- as.data.frame(tmp)
    tmp$d  <- as.numeric(tmp[,nc]) - as.numeric(tmp[,nc-1])
    head(tmp)
    
    g <- tmp %>%
        filter(mc <= 9) %>%
        ggplot(aes(x=time, y=d, colour=factor(mc)))+
        geom_line()+
        facet_wrap(~mc)+
        ggtitle(paste('Difference:',colnames(tmp)[nc],'-', colnames(tmp)[nc-1]))
    plot(g)
}

calc_diff_vs_none <- function(df) {
    # Make sure it's clean:
    df <- as.data.frame(df)
    df <- df[complete.cases(df),]
    
    # Difference vs column "none":
    n <- ncol(df) 
    idx  <- which(grepl('none', names(df)))
    idx2 <- which( (!grepl('none', names(df))) & names(df)!='mc')
    df$diff     <- df[,idx2]-df[,idx]
    df$rel.diff <- (df[,idx2]-df[,idx]) / df[,idx]
    return(df)
}

stats_diff_vs_none <- function(df) {
    return(
        list(mean = mean(df$rel.diff, na.rm=TRUE),
             median = median(df$rel.diff, na.rm=TRUE),
             range = range(df$rel.diff, na.rm=TRUE))
    )
}

comp_dos <- function(z) {
    df <- z %>%
        filter(type=='patient' & wasColonized) %>%
        # DoS per patient:
        group_by(vax.strategy, mc, uid) %>%
        summarize(max.dos = max(DoS, na.rm=TRUE)) %>%
        # Total iso days per MC:
        group_by(vax.strategy, mc) %>%
        summarise(sum.dos = sum(max.dos)) %>%
        # Spread before calculating difference:
        spread(vax.strategy, sum.dos)  %>%
        calc_diff_vs_none()
    
    print('Diff vs none: Duration of stay')
    ss <- stats_diff_vs_none(df)
    print(ss)
    
    title  <-  paste('Duration of stay (rel. diff. vs none)')
    subtitle <- paste('mean:', round(ss$mean,2),
                      '; median:', round(ss$median,2))
    
    g <- ggplot(df)+
        geom_histogram(aes(x=rel.diff), bins=20)+
        geom_vline(xintercept = 0, colour='orange', size=2)+
        ggtitle(title, subtitle)+
        theme(panel.grid = element_blank())
    return(g)
}


plot_histogram_comp <- function(df, 
                                title='', 
                                subtitle='', 
                                x.axis.title='') {
    # DEBUG
    # title='Isolation duration' 
    # subtitle='dummy'
    # x.axis.title='relative difference'
    # # -
    g <- ggplot(df)+
        geom_histogram(aes(x=rel.diff),
                       fill = 'darkgrey',
                       bins=20)+
        geom_vline(xintercept = 0, 
                   colour='black', 
                   size=2,
                   linetype = 'dashed')+
        scale_x_continuous(labels = scales::percent)+
        ggtitle(title, subtitle)+
        xlab(x.axis.title)+
        ylab('Number of simulations')+
        theme(panel.grid = element_blank(),
              text = element_text(size = 22),
              axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 20, l = 0)),
              axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
              plot.margin = unit(x = c(1,1,1,1), units = 'cm'))
    return(g)
}

comp_isolation <- function(z, t.vax.start = 1825) {
    
    # Trim unnecessary data:
    df <- z[z$type=='patient', c('time', 'uid', 'vax.strategy', 'mc', 
                                 'timeIsolationStart','timeIsolationEnd')]
    df  <- df[df$time>t.vax.start]
    
    # Calculate isolation length for each patient:
    df$iso.length <- df$timeIsolationEnd - df$timeIsolationStart
    
    # Retrieve vax scenario name:
    sn = unique(df$vax.strategy)
    scenvaxname = sn[!grepl('none',sn)]
    
    # Summarize by scenario:  
    df <- df %>%
        group_by(vax.strategy, mc) %>%
        summarise(sum.iso = sum(iso.length, na.rm=TRUE)) %>%
        spread(vax.strategy, sum.iso) 
    
    # Identify column of no vax scenario:
    novaxcol <- which(grepl('none',names(df)))
    nc <- ncol(df)
    vaxcol <- ifelse(novaxcol==nc,nc-1,nc)
    zv   = df[,vaxcol] %>% unlist() %>% unname()
    znv  = df[,novaxcol] %>% unlist() %>% unname()
    df$diff.iso = zv - znv  
    df$rel.diff.iso = zv/znv  - 1
    
    # Save result to file:
    write.csv(x = as.data.frame(df), 
                file = paste0('out-compare-iso-',scenvaxname,'.csv'),
                quote = F, row.names = F)
    message(paste('Isolation comparison file saved for',scenvaxname))
    
    g <- df %>%
        ggplot(aes(x=rel.diff.iso)) + 
        geom_histogram(bins=30)+
        ggtitle('Relative change in total isolation days', scenvaxname)
    return(g)
}

comp_sympt <- function(z) {
    
    df <- z %>%
        filter(type == 'patient') %>%
        group_by(vax.strategy, mc, uid) %>%
        summarise(sympt = max(symptomatic, na.rm=TRUE)) %>%
        group_by(vax.strategy, mc) %>%
        summarise(ns = sum(sympt)) %>%
        ungroup() %>%
        spread(vax.strategy, ns) %>%
        calc_diff_vs_none()
    
    print('Diff vs none: symptomatic cases')
    print(df)
    ss <- stats_diff_vs_none(df)
    print(ss)
    
    title  <-  paste('Symptomatic cases (rel. diff. vs none)')
    subtitle <- paste('mean:', round(ss$mean,2),
                      '; median:', round(ss$median,2))
    
    g <- ggplot(df)+
        geom_histogram(aes(x=rel.diff), bins=20)+
        geom_vline(xintercept = 0, colour='orange', size=2)+
        ggtitle(title, subtitle)+
        theme(panel.grid = element_blank())
    return(g)
}

#' Compare symptomatic cases for a given year
comp_sympt_year <- function(z, yr) {
    
    # Add year information: 
    z$year <- 1 + floor(z$time/365)
    
    # Check user error:
    yr.max <- max(z$year)
    if(yr > yr.max){
        msg <- 'ERROR: Cannot compare symptomatic cases for a year that does not exist!'
        msg2 <- paste('; Max year =',yr.max, '; Year requested:',yr)
        stop()
    }
    
    df <- z %>%
        filter(type == 'patient') %>%
        filter(year == yr) %>%
        group_by(vax.strategy, mc, uid) %>%
        summarise(sympt = max(symptomatic, na.rm=TRUE)) %>%
        group_by(vax.strategy, mc) %>%
        summarise(ns = sum(sympt)) %>%
        ungroup() %>%
        spread(vax.strategy, ns) %>%
        calc_diff_vs_none()
    
    print(paste('Diff vs none: symptomatic cases for year',yr))
    print(df)
    ss <- stats_diff_vs_none(df)
    print(ss)
    
    title  <-  paste('Symptomatic cases (rel. diff. vs none) for year',yr)
    subtitle <- paste('mean:', round(ss$mean,2),
                      '; median:', round(ss$median,2))
    
    g <- ggplot(df)+
        geom_histogram(aes(x=rel.diff), bins=20)+
        geom_vline(xintercept = 0, colour='orange', size=2)+
        ggtitle(title, subtitle)+
        theme(panel.grid = element_blank())
    return(g)
}

comp_colectomy <- function(z) {
    df <- z %>%
        filter(type == 'patient') %>%
        group_by(vax.strategy, mc, uid) %>%
        summarise(colec = max(colectomy, na.rm=TRUE)) %>%
        group_by(vax.strategy, mc) %>%
        summarise(nc = sum(colec)) %>%
        ungroup() %>%
        spread(vax.strategy, nc) %>%
        calc_diff_vs_none()
    
    print('Diff vs none: colectomy cases')
    print(df)
    ss <- stats_diff_vs_none(df)
    print(ss)
    
    title  <-  paste('Colectomy cases (rel. diff. vs none)')
    subtitle <- paste('mean:', round(ss$mean,2),
                      '; median:', round(ss$median,2))
    
    g <- ggplot(df)+
        geom_histogram(aes(x=rel.diff), bins=20)+
        geom_vline(xintercept = 0, colour='orange', size=2)+
        ggtitle(title, subtitle)+
        theme(panel.grid = element_blank())
    return(g)
}

ts_prop_vax <- function(z) {
    bck <- 30
    df <- z %>%
        filter(vax.strategy != 'none') %>%
        mutate(tb = round(time/bck)*bck) %>%
        group_by(tb) %>%
        summarise(y = mean(isVax)) %>% 
        mutate(yr = tb/365) 
    
    g <- ggplot(df) +
        geom_line(aes(x=yr, y=y))+
        ggtitle('Proportion of patients vaccinated')+
        xlab('time (year)')+ylab('')
    return(g)
}

#' Proportion of symptomatic cases among vaccinated patients.
ts_prop_sympt_vax <- function(z) {
    bck <- 30
    df <- z %>%
        filter(vax.strategy != 'none') %>%
        filter(isVax) %>%
        mutate(tb = round(time/bck)*bck) %>%
        group_by(tb) %>%
        summarise(y = mean(symptomatic)) %>% 
        mutate(yr = tb/365) 
    
    g <- ggplot(df) +
        geom_line(aes(x=yr, y=y))+
        ggtitle('Proportion of symptomatic among vaccinated')+
        xlab('time (year)')+ylab('')
    g
}

#' Proportion of symptomatic cases among vaccinated patients.
ts_prop_sympt <- function(z) {
    bck <- 30
    df <- z %>%
        filter(wasColonized) %>%
        mutate(tb = round(time/bck)*bck) %>%
        group_by(tb, vax.strategy) %>%
        summarise(y = mean(symptomatic)) %>% 
        mutate(yr = tb/365) 
    
    g <- ggplot(df) +
        geom_line(aes(x=yr, y=y, colour = vax.strategy),
                  size = 2, alpha=0.7)+
        ggtitle('Proportion of symptomatic cases among colonized')+
        xlab('time (year)')+ylab('')
    g
}

sIdx_sympt <- function(z) {
    
    g <- z %>%
        filter(type=='patient') %>%
        ggplot()+
        geom_density(data=z, 
                     aes(x=sIdx, fill= symptomatic),
                     colour = NA,
                     alpha=0.5) +
        facet_wrap(~vax.strategy)
    g
}

OLD_comp_prop_asymptom_dischrg_OLD <- function(z) {
   
    df.pa <- prop_asymptom_dischrg(z)
    
    g <- df.pa %>%
        ggplot() +
        geom_errorbar(aes(x=vax.strategy, 
                          ymin=mn-sd,
                          ymax=mn+sd),
                      size=1) +
        geom_point(aes(x=vax.strategy, y=mn), size=4)+
        geom_text(aes(x=vax.strategy, y=mn, label=round(mn,5)),nudge_x = 0.1)+
        ggtitle('Proportion asymptomatic at discharge time',
                'Mean +/- sd')+
        ylab('Proportion')
    g
}

#' Proportion of asymptomatic cases at discharge time
comp_prop_asymptom_dischrg <- function(z, t.vax.start = 1825) {
    
    # Scenario name:
    sn <- unique(z$vax.strategy)
    scenvaxname <- sn[!grepl('none',sn)]
    
    message(paste('Comparing asymptomatic cases for',scenvaxname,'...'))
    
    # Among colonized only: 
    z <- z[z$wasColonized,]
    
    # Only necessary data:
    df <- z[z$time>t.vax.start, c('mc', 'symptomatic','vax.strategy')]
    
    # Calculate grouped numbers & proportions:
    dfa <- df %>%
        group_by(mc,vax.strategy) %>%
        summarise(n.asympt = sum(!symptomatic),
                  p.asympt = sum(!symptomatic) / n())
    
    # Calculate differences for numbers:
    dfa.n <- dfa %>%
        select(-p.asympt) %>%
        spread(vax.strategy, n.asympt)
    
    col.vax   <- which(grepl('frailty',names(dfa.n)))
    col.novax <- which(grepl('none',names(dfa.n)))
    
    v0 <- dfa.n[,col.novax] %>% unlist(use.names = F)
    v1 <- dfa.n[,col.vax] %>% unlist(use.names = F)
    
    dfa.n$diff <- v1-v0
    dfa.n$rel.diff <- v1/v0 - 1
    dfa.n
    
    # Calculate differences for numbers:
    dfa.p <- dfa %>%
        select(-n.asympt) %>%
        spread(vax.strategy, p.asympt)
    
    col.vax   <- which(grepl('frailty',names(dfa.p)))
    col.novax <- which(grepl('none',names(dfa.p)))
    
    v0 <- dfa.p[,col.novax] %>% unlist(use.names = F)
    v1 <- dfa.p[,col.vax] %>% unlist(use.names = F)
    
    dfa.p$diff <- v1-v0
    dfa.p$rel.diff <- v1/v0 - 1
    dfa.p
    
    # save CSV files:
    fname.n <- paste0('out-compare-asympt-num-among-colonized-',scenvaxname,'.csv')
    write.csv(x = dfa.n, file = fname.n, quote = F, row.names = F)
    fname.p <- paste0('out-compare-asympt-prop-among-colonized-',scenvaxname,'.csv')
    write.csv(x = dfa.p, file = fname.p, quote = F, row.names = F)
    
    message(paste("Asymptomatic cases among colonized: CSV file saved for",scenvaxname))
    
    # Plot    
    g <- dfa.p %>%
        ggplot(aes(x=rel.diff)) + 
        geom_histogram(bins=30)+
        ggtitle("Change in prop. of asymptomatic cases among colonized", scenvaxname)
    return(g)
}


comp_environment_acq <- function(z, t.vax.start=1825) {
    
    sn <- unique(z$vax.strategy)
    vaxscen <- sn[grepl('frailty', sn)]

    # After vax:
    z <- z[z$time > t.vax.start,]
    
    # filter acquisition that occured in hospital only:
    a <- z[z$acquisitionType=='environment' | z$acquisitionType=='hcw_asymptomatic',]
    
    # counts:
    a <- a %>%
        group_by(mc, vax.strategy, acquisitionType) %>%
        tally() %>%
        spread(acquisitionType, n)
    
    # Proportion of environmental acquisition:
    a$prop.env <- a$environment / (a$environment+a$hcw_asymptomatic)
        
    # Save CSV file:
    write.csv(x = a, 
              file = paste0('out-colonized-hai-prop-env-',
                            vaxscen,
                            '.csv'),
              quote = F, row.names = F)
    
}

comp_acqType <- function(z) {
    a <- z %>%
        group_by(vax.strategy,mc, uid, acquisitionType) %>%
        summarise(n = n_distinct(acquisitionType)) %>%
        filter(n==1) %>%
        group_by(vax.strategy,acquisitionType) %>%
        summarise(nat = sum(n)) %>%
        filter(acquisitionType != 'NA')%>%
        mutate(p = nat/sum(nat)) 

    g <- a %>%    
        ungroup() %>%
        ggplot(aes(x=acquisitionType, y=p, fill=vax.strategy)) +
        geom_bar(stat='identity', position = 'dodge')
    g
}

comp_prop_symptom <- function(z) {
    summary(z$time)
    
    a <- z %>%
        filter(wasColonized) %>%
        filter(time>2400) %>%
        group_by(vax.strategy, uid,mc) %>%
        summarise(x = max(symptomatic)) %>%
        group_by(vax.strategy) %>%
        summarise(p = mean(x))
    print(a)
    res <- a$p[a$vax.strategy=='frailty'] / a$p[a$vax.strategy=='none']
    return(res)
}

#' Number of vaccines required to prevent one symptomatic case
vax_req_prevent_one_case <- function(z, t.start.vax = 1825){
    
    # Scenario name:
    sn <- unique(z$vax.strategy)
    scenvaxname <- sn[!grepl('none',sn)]
    
    message("Number of vax required to prevent one symptomatic case for:")
    message(scenvaxname)
    
    # Filter necessary variables only:
    df <- z[z$time> t.start.vax, 
            c('mc', 'isVax', 'symptomatic', 'vax.strategy')]
    
    # Summarize data:
    df <- df %>%
        group_by(vax.strategy, mc) %>%
        summarize(n.s = sum(symptomatic), 
                  n.vax = sum(isVax))
    
    # Calculate difference of number of symptomatic:
    df1 <- select(df, -n.vax) %>%
        spread(vax.strategy, n.s) 
    
    k = which(names(df1)==scenvaxname)
    nc = ncol(df1)
    k0 = ifelse(k==nc, nc-1,nc)
    
    df1$diff.n.s <- unlist(df1[,k0]) - unlist(df1[,k])
    
    # Calculate number of vaccines required 
    # to prevent one symptomatic case:
    df2 <- df %>% 
        filter(vax.strategy==scenvaxname) %>%
        left_join(df1, by='mc') %>%
        mutate(vax.prev.one = n.vax/diff.n.s)
    
    # Save to CSV file:
    fname = paste0('out-vax-prev-one-',scenvaxname,'.csv')
    write.csv(file = fname, x = df2, quote = F, row.names = F)
    
    message(paste("CSV file saved:", fname))
    
    # plot 
    g <- df2 %>%
        ggplot(aes(x=vax.prev.one)) + 
        geom_histogram(bins=30) + 
        ggtitle("Number of vax required to prevent one symptomatic case", 
                scenvaxname)
    return(g)
}


# ---- Main Function -----

full_comparison <- function(scen.pair) {   # scen.pair = c('none-0-prm_vax', 'frailty-80-prm_vax_p_0')
    
    scenvax <- scen.pair[!grepl('none',scen.pair)]
    
    print(paste('Loading all simulations RData for',
                scenvax,'...'))
    z <- load_all_simul(scen.pair)
    print('=> All simulations RData loaded. ')
    
    # Compare variables:
    print(" Start comparing variables...")
    g.iso   <- comp_isolation(z)
    g.dos   <- comp_dos(z)
    g.asymp <- comp_prop_asymptom_dischrg(z)
    g.vaxprv<- vax_req_prevent_one_case(z, t.start.vax = 1825)
    g.env   <- comp_environment_acq(z)
    print("=> all comparisons processed.")
    
    # Plots
    print('Start plotting comparisons...')
    fname = paste0('plot-compare-',scenvax,'.pdf')
    pdf(file = fname, 
        width = 18, height = 12)
    grid.arrange(g.asymp, g.vaxprv,
                 g.iso, g.dos,
                 ncol=2)
    dev.off()
    print('=> plots completed.')
}


# ---- Run ----

t1 <- as.numeric(Sys.time())

# Select strategies:
scen.pair <- list(
    c('none-0-prm_vax', 'frailty-40-prm_vax_p_0'),
    c('none-0-prm_vax', 'frailty-80-prm_vax_p_0')
)


message(" >>> Comparing scenarios ...")
message(paste(scen.pair, collapse = ' vs. '))

# Run the comparisons:
x <- lapply(scen.pair, full_comparison)

# End stuff
t2 <- as.numeric(Sys.time())
dt <- t2-t1
dtm <- round(dt/60,1)
msg <- paste0('== Comparison script done in ',dtm,' min ==')
message(msg)
print(msg)


# --- old stuff:
# g.p.vax <- ts_prop_vax(z)
# g.p.svax<- ts_prop_sympt_vax(z)
# g.psympt<- ts_prop_sympt(z)
# g.sidx  <- sIdx_sympt(z)
# g.pad   <- comp_prop_asymptom_dischrg(z)
# 
# ---