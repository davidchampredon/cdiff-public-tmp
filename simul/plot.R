###
###  CHECKS ON SIMULATION RUN
###

library(tidyverse)

# Load simulation output objects:
load('simul.RData')

pdf('plot.pdf', width = 15, height = 10)

theme_set(theme_bw())

# ---- Susceptibility ----
df %>%
    group_by(uid, type, hcs_name) %>%
    summarise(s = mean(sIdx)) %>%
    ggplot() +
    geom_histogram(aes(x=s, fill=type), alpha=0.8, binwidth = 0.05)+
    ggtitle('Susceptibility Index') +
    coord_cartesian(xlim=c(0,1)) + guides(fill=FALSE)+
    facet_grid(hcs_name~type, scales = 'free')+
    xlab('susceptibility index')

df %>%
    filter(isVax) %>%
    filter(uid <25)  %>%
    ggplot(aes(x=time, y=sIdx)) +
    geom_line(size=2)+
    facet_wrap(~uid, scales = 'fixed')+
    ggtitle('Susceptibility index','selected vaccinated patients')+
    ylab('Susceptibility index')


# ---- DoS ----
df %>%
    filter(type=='patient') %>%
    group_by(uid) %>%
    summarize(d = max(DoS)) %>%
    ggplot() +
    geom_histogram(aes(x=d),fill='darkgrey', bins=30) +
    scale_x_log10(breaks=c(1,2,5,7,10,30,90,365,730))+
    ggtitle("Duration of stay (days)")

# ---- Incubation ----
df %>%
    filter(symptomatic==TRUE) %>%
    group_by(uid) %>%
    summarize(d = mean(prdIncubation)) %>%
    ggplot() +
    geom_histogram(aes(x=d),fill='blue', bins = 12) +
    ggtitle("Incubation period")

# ---- Symptomatic period ----
df %>%
    filter(symptomatic==TRUE) %>%
    group_by(uid) %>%
    summarize(d = mean(prdSymptom)) %>%
    ggplot() +
    geom_histogram(aes(x=d),fill='orange', bins=20) +
    ggtitle("Symptomatic period")

# ---- Population ----

df %>%
    group_by(time, type) %>%
    summarise(n = length(type)) %>%
    ggplot() +
    geom_step(aes(x=time,y=n, colour=type), size=2)+
    ggtitle('Population')

g.adm <- df %>%
    filter(type=='patient') %>%
    group_by(uid) %>%
    summarize(ta = min(timeAdmission)) %>%
    group_by(ta) %>%
    summarize(n = length(uid)) %>%
    filter(ta>0) %>%
    mutate(cum.adm = cumsum(n)) %>% 
    ggplot(aes(x=ta,y=cum.adm)) + 
    geom_step() + 
    geom_smooth(method = 'lm', size=0.2, se = F) + 
    geom_abline(intercept = 0, 
                slope = prm.hosp$admission_rate, 
                linetype=2, colour='red')+
    xlab('time')+
    ggtitle('Cumulative admissions')
plot(g.adm)


# ---- Contact durations ----

x <- res[['contact_duration_Hospital']] * 1440  # 1440 minutes in one day
mc <- mean(x)
hist(x, breaks = 60,
     col = 'lightgrey', border='grey',
     main = paste('HCW/patient contact durations\n mean =',
                  round(mc,1)),
     xlab = 'minutes')
abline(v=mc, lty=2,lwd=2)

xh <- res[['contact_duration_Hospital']] * 1440  # 1440 minutes in one day
xl <- res[['contact_duration_LTCF']] * 1440  # 1440 minutes in one day

df.contact <- data.frame(d = c(xh,xl),
                        hcs = c(rep('Hospital',length(xh)),
                                rep('LTCF',length(xl))))

df.contact %>%
    ggplot()+
    geom_histogram(aes(x=d, fill=hcs), 
                   position='identity',bins=50, alpha=0.5) +
    geom_vline(data=df.contact %>%
                   group_by(hcs) %>%
                   summarize(m=mean(d)),
               aes(xintercept=m, colour=hcs), linetype=2)+
    xlab('duration (minute)')+
    ggtitle('Contact durations HCW/patients')


# ---- Time series of prevalence ----
ts.c <- df %>%
    group_by(time, full_type) %>%
    summarise(n = sum(isColonized))
ggplot(ts.c, aes(x=time, y=n, colour=full_type)) +
    geom_step(size=3, alpha=0.6) + #geom_point()+
    ggtitle('Colonized')

# by room
ts.c.r <- df %>%
    filter(type == 'patient') %>%
    group_by(time, hcs_name, current_room_uid) %>%
    summarise(n.patients = length(uid),
              n.colonized = sum(isColonized),
              n.symptomatic = sum(symptomatic)) %>%
    gather('type','value',4:6)

g.ts <- ggplot(ts.c.r, aes(x=time, y=value, colour=type)) +
    geom_point(alpha=0.7,size=2) +
    facet_wrap(~hcs_name+current_room_uid, 
               scales = 'free_y')+
    ggtitle('Patient occupancy, Colonized & Symptomatic per room')
    

# ---- Contamination ----

# Rooms:
rc <- df %>%
    filter(current_room_uid < 9999) %>%
    group_by(time, hcs_name, current_room_uid) %>%
    summarise(el = mean(current_room_contamIdx),
              n = length(current_room_contamIdx))
ggplot(rc, aes(x=time, y=el)) +
    geom_line(size=0.5, colour='gray') + geom_point(size=0.5) +
    facet_wrap(~hcs_name+current_room_uid) +
    ggtitle('Room environmental contamination') + ylab('Load')


# Individuals' contamination index:
# (only for a sample of the population, else too large)

icitry <- try(
    ici <- df %>%
        filter(isColonized) %>%
        filter(uid %in% sample(unique(uid),20,replace = FALSE)) %>%
        select(hcs_name, full_type, uid, time, contamIdx, currentlySymptomatic) %>%
        ggplot(aes(x=time, y=contamIdx, colour=currentlySymptomatic)) + 
        geom_point()+ geom_line(alpha=0.4) +
        facet_wrap(hcs_name ~ full_type + uid, scales='free')+
        ggtitle('Individuals contamination index'),
    silent = TRUE)

if(class(icitry)[1] == 'try-error') 
    plot(ggplot()+ggtitle('UNSUCCESSFUL PLOT (sampling colonized individuals)')) 
if(class(icitry)[1] != 'try-error') 
    plot(ici)
    
    
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
    ggtitle('Other environmental contamination (dashed: Mean contamination index of all patient\'s rooms)')

# ---- Shedding index ----
if(FALSE){
ts.shi <- df %>%
    group_by(time, uid, type) %>%
    summarise(shi = mean(contamIdx)) %>%
    filter(uid %in% sample(unique(uid),20,replace = FALSE)) %>%
ggplot(ts.shi,aes(x=time, y=shi, colour=factor(type))) +
    geom_line() + #geom_point()+
    ggtitle('contamIdx') +
    facet_wrap(~uid, scales='free_y') +
    guides(colour=FALSE)
}

# ---- Vaccination ----

vax <- df %>%
    group_by(hcs_name,type) %>%
    summarise(pvax = sum(isVax)/length(isVax))

df %>%
    group_by(hcs_name,type,time) %>%
    summarise(nvax = sum(isVax)) %>%
    ggplot(aes(x=time, y=nvax)) +
    geom_line(aes(colour=type), size=1)+
    facet_wrap(~hcs_name)+
    ggtitle('number of vaccinated individuals')




# ---- End ----
dev.off()
message(' - ALL PLOTS DONE - ')


