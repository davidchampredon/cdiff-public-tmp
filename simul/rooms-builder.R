library(tidyr); library(dplyr); library(ggplot2)
#library(tidyverse)

set.seed(1234)

#' Create CSV file of rooms descriptions
#' @param n.iso Integer. Number of isolation rooms.
#' @param n.wards Integer. Number of wards.
#' @param rooms.per.ward.mean Numeric. Mean number of rooms per ward (Poisson distribution).
#' @param beds.per.room.mean Numeric. Mean number of beds per room (Poisson distribution).
#' @param filename String. Name of the CSV file to save the room data to. 
build_rooms <- function(n.iso, 
                        n.wards, 
                        rooms.per.ward.mean, 
                        beds.per.room.mean,
                        min.bed.per.room,
                        max.bed.per.room,
                        filename) {
    
    # Isolation rooms:
    df.iso <- data.frame(type='isolation_room',
                         ward = 0, 
                         name = paste0('ISO',1:n.iso),
                         max_patient = 1)
    
    n.rooms.per.ward <- rpois(n=n.wards, lambda = rooms.per.ward.mean)
    n.rooms.per.ward[n.rooms.per.ward==0] <- 1
    
    # Patient rooms
    df.pr <- list()
    
    for(w in 1:n.wards){
        
        lmbd  <- beds.per.room.mean - min.bed.per.room
        mbeds <- rpois(n = n.rooms.per.ward[w], lambda = lmbd)
        
        n.beds.per.room <- min.bed.per.room + mbeds
        
        # apply limits:
        n.beds.per.room[n.beds.per.room < min.bed.per.room] <- min.bed.per.room
        n.beds.per.room[n.beds.per.room > max.bed.per.room] <- max.bed.per.room
        
        df.pr[[w]] <- data.frame(type='patient_room',
                                 ward = w, 
                                 name = paste0(w,'_',LETTERS[1:n.rooms.per.ward[w]]),
                                 max_patient = n.beds.per.room)
    }
    df.pat.room <- do.call('rbind', df.pr)
    
    # Binds isolation and patient rooms:
    df.rooms <- rbind(df.iso, df.pat.room)
    
    # Save to CSV file:
    write.csv(df.rooms, file = filename,
              row.names = FALSE, quote = FALSE)
    
    # Return object
    return(list(df = df.rooms, 
                filename = filename))
}

summary_rooms <- function(x) {
 
    n <- sum(x$df$max_patient)
    
    print(paste0('There is a maximum of ',n, ' patients in ',x$filename))
    print('Maximum patients per ward:')
    x$df %>%
        group_by(ward) %>%
        summarise(n = sum(max_patient)) %>%
        print()   
}

plot_rooms <- function(x){
    
    df <- x$df
    
    g <- df %>%
        mutate(is.ISO = ifelse(grepl('ISO',name),'Isolation room','Standard room')) %>%
        ggplot()+
        geom_histogram(aes(x=max_patient), 
                       fill = 'grey',
                       colour='darkgrey',
                       binwidth=1) +
        facet_wrap(~is.ISO)+
        theme(panel.grid = element_blank())+
        ggtitle('Room sizes (overall)')
    
    gw <- df %>%
        mutate(is.ISO = grepl('ISO',name)) %>%
        filter(!is.ISO) %>%
        ggplot()+
        geom_histogram(aes(x=max_patient), 
                       fill = 'grey',
                       colour='darkgrey',
                       binwidth=1) +
        facet_wrap(~ward, scales = 'fixed')+
        theme(panel.grid = element_blank())+
        ggtitle('Room size by ward')
    
    pdf('plot-rooms.pdf')
    plot(g)
    plot(gw)
    dev.off()
}



do.test <- 0

if(do.test){
    
    r.hosp <- build_rooms(n.iso   = 80 , 
                          n.wards = 10 ,
                          rooms.per.ward.mean = 11, 
                          beds.per.room.mean  = 2.1,
                          min.bed.per.room = 2,
                          max.bed.per.room = 4,
                          filename = 'rooms_test.csv')
   
    summary_rooms(r.hosp)
}



