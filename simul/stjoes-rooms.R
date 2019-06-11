library(tidyverse)
dat <- read.csv('stjoes-rooms.csv')


n  <- sum(dat$census)
ns <- sum(dat$singleRoom)
ns/n

# Anne Bialachowski: 
# The vast majority of rooms that have more than 1 bed have 2. 
# There are probably only 2, 4 bed rooms in the building.
mean.n.beds.per.room <- 2.1


dat2 <- dat %>%
    filter(!grepl('ICU',Unit)) %>%
    mutate(patients.in.multibed = census - singleRoom) %>%
    mutate(n.rooms = patients.in.multibed / mean.n.beds.per.room)

mean(dat2$n.rooms)

# hist(dat2$multibed)
