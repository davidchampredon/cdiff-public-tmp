###
###  CHECKS ON SIMULATION RUN
###

library(tidyverse)
theme_set(theme_bw())

# Load simulation output objects:
load('simul-mc.RData')

# Only deal with one HCS:
stopifnot(length(unique(df$hcs_name))==1)



# ---- Admissions ----

df.adm <- df %>%
    group_by(timeAdmission,mc)%>%
    summarise(n = n_distinct(uid)) %>%
    mutate(yr = 1+timeAdmission%/%365)

df.adm %>%
    group_by(yr,mc) %>%
    summarise(adm.yr = sum(n)) %>%
    ggplot(aes(x=factor(yr), y=adm.yr)) +
    geom_boxplot(fill='grey')+
    ggtitle('Annual number of admissions')+
    xlab('year') + ylab('number of admissions')


# ---- Patient Days ----

dt <- max(diff(df$time))

pd <- df %>% 
    filter(type=='patient') %>%
    mutate(month = 1 + time%/%30) %>%
    group_by(wardAssigned, month,mc) %>%
    summarise(patient.days = round(n()*dt))

pd %>%
    ggplot(aes(x=month, y=patient.days, colour=factor(mc)))+
    geom_point()+geom_line()+
    facet_wrap(~wardAssigned)+
    ggtitle('Monthly patient days, by ward')

pd %>%
    group_by(month, mc) %>%
    summarise(sum.pd = sum(patient.days)) %>%
    group_by(month) %>%
    summarise(mean.pd = mean(sum.pd),
              sd.pd = sd(sum.pd)) %>%
    ggplot(aes(x=month, y=mean.pd))+
    geom_line()+
    geom_ribbon(aes(ymin=mean.pd-sd.pd, ymax=mean.pd+sd.pd),
                alpha=.2)+
    ggtitle('Mean (sd) monthly patient-days for the whole HCS (across all MCs)')
