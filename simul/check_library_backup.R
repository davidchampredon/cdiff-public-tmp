library(dplyr)
library(tidyr)
library(ggplot2)

R.library.dir    <- '../Rlibrary/lib'
library(abmcdiff,lib.loc = R.library.dir)


rngseed <- 123456

set.seed(rngseed)
tictac <- as.numeric(Sys.time())

# ---- Set up ----

create.list <- function(z) {
    param <- list()
    for(i in 1:nrow(z)){
        if(is.na(as.numeric(z$value[i])))
            param[[z$name[i]]] <- z$value[i]
        else{
            param[[z$name[i]]] <- as.numeric(z$value[i])
        }
    }
    return(param)
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

prm <- full.param.HCS(file_infrastruct = 'rooms_hosp.csv', 
                      file_epi         = 'prm_hosp.csv') 

sim.prm <- list(horizon = 200, 
               timestep = 0.5,
               seedRNG  = rngseed) # IMPORTANT : do not forget "set.seed()" in this R script

disease.prm <- read.csv('prm_disease.csv', 
                        as.is = T, strip.white = T) %>%
    create.list()


# --- Run Simulations ----

res <- abmcdiff_test_4(prm, sim.prm, disease.prm)

res.pop <- res[grepl('pop_',names(res))]
df.pop.dt <- lapply(X = res.pop, FUN = as.data.frame)
df <- do.call(what = 'rbind', args = df.pop.dt)

df$full_type <- paste(df$type, df$subtype, sep='_')

df.ward <- df %>% 
    group_by(time,wardAssigned) %>% 
    summarize(n=length(uid))

# --- Plots ----

pdf('plot.pdf', width = 15, height = 10)

df %>% 
    group_by(uid, type) %>%
    summarise(s = mean(sIdx)) %>%
    ggplot() + 
    geom_histogram(aes(x=s, fill=type), alpha=0.8, binwidth = 0.05)+ 
    ggtitle('Susceptibility Index') + 
    coord_cartesian(xlim=c(0,1)) + guides(fill=FALSE)+
    facet_wrap(~type, scales = 'free')

# DoS
df %>% 
    filter(type=='patient') %>% 
    group_by(uid) %>%
    summarize(d = max(DoS)) %>%
    ggplot() + 
    geom_histogram(aes(x=d),fill='darkgrey', binwidth = 1) + 
    ggtitle("Duration of stay")

# Incubation
df %>% 
    filter(symptomatic==TRUE) %>% 
    group_by(uid) %>%
    summarize(d = mean(prdIncubation)) %>%
    ggplot() + 
    geom_histogram(aes(x=d),fill='blue', binwidth = 1) + 
    ggtitle("Incubation period")

# Symptomatic period
df %>% 
    filter(symptomatic==TRUE) %>% 
    group_by(uid) %>%
    summarize(d = mean(prdSymptom)) %>%
    ggplot() + 
    geom_histogram(aes(x=d),fill='orange', binwidth = 1) + 
    ggtitle("Symptomatic period")


# Population
dfs <- df %>% 
    group_by(time) %>%
    summarise(np = sum(type=='patient'),
              nh = sum(type=='HCW')) 
ggplot(dfs,aes(x=time, y=np)) + geom_step() + geom_point() + geom_point(aes(y=nh), color='green')



# Time series of prevalence
ts.c <- df %>% 
    group_by(time, full_type) %>%
    summarise(n = sum(isColonized))
ggplot(ts.c, aes(x=time, y=n, colour=full_type)) + 
    geom_step(size=3, alpha=0.6) + #geom_point()+
    ggtitle('Colonized')

# by room
ts.c.r <- df %>% 
    filter(type == 'patient') %>%
    group_by(time, current_room_uid) %>%
    summarise(n = length(uid),
              nc = sum(isColonized),
              ns = sum(symptomatic))
ggplot(ts.c.r, aes(x=time)) + 
    geom_step(aes(y=n),size=0.5, colour='black') +
    geom_step(aes(y=nc),size=1.2, colour='orange') +
    geom_step(aes(y=ns),size=1.2, colour='red') +
    facet_wrap(~current_room_uid, scales = 'free_y') + 
    ggtitle('Patient occupancy, Colonized & Symptomatic per room')

# Room contamination:
rc <- df %>% 
    filter(current_room_uid < 9999) %>%
    group_by(time, current_room_uid) %>%
    summarise(el = mean(current_room_contamIdx), 
              n = length(current_room_contamIdx))
ggplot(rc, aes(x=time, y=el)) + 
    geom_line(size=0.5, colour='gray') + geom_point(size=0.5) + 
    facet_wrap(~current_room_uid) + 
    ggtitle('Room environmental contamination') + ylab('Load')


# individual contamination index:

if(length(unique(df$uid))<200){
    
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

# Other environmental contamination:
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
              
              

# Shedding index
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

# Check "on duty" schedules:
duty <- df %>%
    filter(type=='HCW') %>%
    filter(time < 20) %>%
    group_by(uid, time) %>%
    summarise(y = mean(onDuty))

ggplot(duty, aes(x=time, y=as.logical(y))) +
    geom_point() +
    facet_wrap(~uid) +
    ggtitle('HCW schedule (20 first days)') + ylab('On Duty')


# Contact durations
x <- res[['contact_duration']] * 1440  # 1440 minutes in one day
mc <- mean(x)
hist(x, breaks = 60,
     col = 'lightgrey', border='grey',
     main = paste('HCW/patient contact durations\n mean =',
                  round(mc,1)),
     xlab = 'minutes')
abline(v=mc, lty=2,lwd=2)

dev.off()

tictac2 <- as.numeric(Sys.time())
dt <- tictac2 - tictac
message(paste(" - - - DONE IN ", dt%/%60,'min ', round(dt%%60,0), "sec - - - "))
message(' - - - ')